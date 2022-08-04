import os
import rasterio
import rioxarray
import geopandas as gpd
import subprocess
import fiona
import fiona.crs
import pyvista as pv
import shapely.geometry
import shapely.geometry.polygon
from tin_maker.utils import get_output_dir, polygon_to_ring, extract_segments, \
    triangle2shp, get_largest_shape


class TinGenerator(object):

    def __init__(self, watershed_tif, watershed_shp, streams_shp, outlet_csv):

        ds = rasterio.open(watershed_tif)
        self.ds = ds

        affine = ds.meta['transform']
        self.cell_x, self.cell_y = affine[0], -1 * affine[4]

        rds = rioxarray.open_rasterio(watershed_tif)
        rds = rds.squeeze().drop("spatial_ref").drop("band")
        rds.name = "data"

        df = rds.to_dataframe().reset_index()
        df['geometry'] = df.apply(lambda x: shapely.geometry.Point(x['x'], x['y']), axis=1)

        sub_df = df[df.data >= 0.0]
        gdf = gpd.GeoDataFrame(sub_df, geometry=sub_df['geometry'], crs=ds.crs)

        self.gdf = gdf
        self.watershed_boundary = get_largest_shape(watershed_shp)
        self.streams_gdf = gpd.read_file(streams_shp)

        # set outlet
        self.outlet_x = -1
        self.outlet_y = -1
        with open(outlet_csv, 'r') as fptr:
            outlet_line = fptr.read()
            outlet_line_parts = outlet_line.split(",")
            self.outlet_x = float(outlet_line_parts[0])
            self.outlet_y = float(outlet_line_parts[1])

    def _decimate_nodes(self, decimation_factor=0.25,
                       watershed_boundary_buffer_dist=0.0,
                       stream_boundary_buffer_dist=0.0):

        if watershed_boundary_buffer_dist > 0:
            watershed_boundary_ring = polygon_to_ring(self.watershed_boundary, watershed_boundary_buffer_dist)
            watershed_boundary_points = self.gdf.within(watershed_boundary_ring)
            self.gdf = self.gdf.loc[watershed_boundary_points == False]

        if stream_boundary_buffer_dist > 0:
            gdf_streams_buffer = self.streams_gdf.buffer(stream_boundary_buffer_dist, resolution=64)
            dissolved_buffer = gdf_streams_buffer.geometry.unary_union
            stream_points = self.gdf.geometry.within(dissolved_buffer)
            self.gdf = self.gdf.loc[stream_points == False]

        # clean mesh and decimate points (decimation done using pyvista mesh since triangle doesn't decimate by elevation)
        xyz_points = list(zip(self.gdf['x'], self.gdf['y'], self.gdf['data']))
        mesh = pv.PolyData(xyz_points)
        mesh = mesh.clean(tolerance=0.5, absolute=True) #merge points < 0.5 m apart
        surf = mesh.delaunay_2d()
        surf = surf.decimate(decimation_factor)

        # recreate xyz points from decimated mesh
        xyz_points = []
        for pnt in surf.points:
            xyz_points.append((pnt[0], pnt[1], pnt[2]))

        return xyz_points

    def generate_tin(self, odir, decimation_factor,
                       watershed_boundary_buffer_dist=0.0,
                       stream_boundary_buffer_dist=0.0):

        # interior nodes (decimated)
        xyz_points = self._decimate_nodes(decimation_factor, watershed_boundary_buffer_dist, stream_boundary_buffer_dist)

        # channel segments
        channel_segments, channel_points = extract_segments(self.streams_gdf, 'chan', odir)
        channel_points_z = {}

        for i, seg in enumerate(channel_segments):
            x1, y1 = channel_points[seg[0]]
            x2, y2 = channel_points[seg[1]]
            zs = self.ds.sample([(x1, y1), (x2, y2)])
            z1 = next(zs)
            z2 = next(zs)
            channel_points_z[seg[0]] = z1[0]
            channel_points_z[seg[1]] = z2[0]

        for i, pnt in enumerate(channel_points):
            channel_points[i] = pnt[0], pnt[1], channel_points_z[i]

        # watershed boundary segments
        wshed_points_z = {}
        wshed_segments, wshed_points = extract_segments(self.watershed_boundary, 'catch', odir)
        for i, seg in enumerate(wshed_segments):
            x1, y1 = wshed_points[seg[0]]
            x2, y2 = wshed_points[seg[1]]
            zs = self.ds.sample([(x1, y1), (x2, y2)])
            z1 = next(zs)
            z2 = next(zs)
            wshed_points_z[seg[0]] = z1[0]
            wshed_points_z[seg[1]] = z2[0]

        for i, pnt in enumerate(wshed_points):
            wshed_points[i] = pnt[0], pnt[1], wshed_points_z[i]

        nodes_file = os.path.join(odir, 'out.node')
        edges_file = os.path.join(odir, 'out.edge')
        if os.path.exists(nodes_file):
            os.remove(nodes_file)

        fptr1 = open(nodes_file, 'w')
        points = xyz_points + channel_points + wshed_points
        point_types = [0] * len(xyz_points) + [2] * len(channel_points) + [1] * len(wshed_points)

        # set the channel node closest to outlet_x, outlet_y as the outlet node
        outlet_dist = 1e12
        outlet_node_idx = None
        for i, pnt in enumerate(points):
            dist = (pnt[0] - self.outlet_x) * (pnt[0] - self.outlet_x) + \
                   (pnt[1] - self.outlet_y) * (pnt[1] - self.outlet_y)
            if dist < outlet_dist:
                outlet_node_idx = i

        point_types[outlet_node_idx] = 3

        num_vertices = len(points)

        # First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
        fptr1.write(f'{num_vertices} 2 1 1\n')
        for i, v in enumerate(points):
            fptr1.write(f'{i} {v[0]} {v[1]} {v[2]} {point_types[i]}\n')
        fptr1.close()

        # First line: <# of edges> <# of boundary markers (0 or 1)>
        # Following lines: <edge #> <endpoint> <endpoint> [boundary marker]
        fptr2 = open(edges_file, 'w')
        num_edges = len(channel_segments + wshed_segments)
        fptr2.write(f'{num_edges} 1\n')

        num_nodes = len(xyz_points)
        for i, seg in enumerate(channel_segments):
            fptr2.write(f'{i} {seg[0] + num_nodes} {seg[1] + num_nodes} 2\n')

        num_nodes = len(xyz_points) + len(channel_points)
        for i, seg in enumerate(wshed_segments, len(channel_segments)):
            fptr2.write(f'{i} {seg[0] + num_nodes} {seg[1] + num_nodes} 1\n')

        fptr2.close()

        poly_file = os.path.join(odir, 'out.poly')
        if os.path.exists(poly_file):
            os.remove(poly_file)
        with open(poly_file, 'w') as fptr3:
            with open(nodes_file, 'r') as fptr1:
                nodes_data = fptr1.read()
                fptr3.write(nodes_data)

            with open(edges_file, 'r') as fptr2:
                edges_data = fptr2.read()
                fptr3.write(edges_data)

            fptr3.write('0\n')

        return nodes_file, edges_file, poly_file
