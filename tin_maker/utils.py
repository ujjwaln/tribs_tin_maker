import os
import shapely.geometry
from geojson import FeatureCollection, LineString
import geopandas as gpd
import pandas as pd
import json
import fiona
import subprocess
import json
from pysheds.grid import Grid


class Watershed(object):
    def __init__(self, watershed_tif=None, watershed_shp=None, channel_shp=None, outlet_csv=None, area=None, id=None):
        self.watershed_tif = watershed_tif
        self.watershed_shp = watershed_shp
        self.channel_shp = channel_shp
        self.outlet_csv = outlet_csv
        self.area = area
        self.id = id
        if watershed_tif:
            self.folder = os.path.dirname(watershed_tif)

    def from_dict(self, dic):
        self.__dict__.update(dic)
        if self.watershed_tif:
            self.folder = os.path.dirname(self.watershed_tif)

    def toJson(self):
        return json.dumps(self.__dict__)

    def __repr__(self):
        return self.toJson()


def run_triangle(poly_file):
    triangle_path = './triangle/triangle'
    params = '-p'
    result = subprocess.run(
        [triangle_path, params, poly_file], capture_output=True
    )
    # print("stdout:", result.stdout)
    # print("stderr:", result.stderr)



def get_output_dir(name, base_dir='./'):
    base_dir = os.path.abspath(base_dir)
    odir = os.path.join(base_dir, name)
    if not os.path.exists(odir):
        os.makedirs(odir)
    return odir


def geojson2shapely(line):
    assert (isinstance(line, LineString))
    return shapely.geometry.LineString(line["coordinates"])


def polygon_to_ring(boundary: shapely.geometry.Polygon, buffer_distance):
    buffered = boundary.buffer(distance= -1 * buffer_distance)
    # return buffered
    return boundary.difference(buffered)


def get_largest_shape(shapefile):
    shapes = fiona.open(shapefile)
    largest_shape = None
    for shp_dict in shapes:
        shp = shapely.geometry.shape(shp_dict['geometry'])
        if largest_shape is None:
            largest_shape = shp
        elif shp.area > largest_shape.area:
            largest_shape = shp
        else:
            pass
    return largest_shape


def extract_segments(geoms, output_file_prefix=None, output_dir=None):
    coords = []
    segments = []
    if isinstance(geoms, str) or isinstance(geoms, gpd.GeoDataFrame):
        if isinstance(geoms, str):
            geoms_gdf = gpd.read_file(geoms)
        else:
            geoms_gdf = geoms

        for index, row in geoms_gdf.iterrows():
            geom = row.geometry
            if geom.geom_type == 'Polygon':
                j = len(coords)
                coords += geom.exterior.coords
                for i in range(j, j+len(geom.exterior.coords)-1):
                    segments.append((i, i+1))
            elif geom.geom_type == 'LineString':
                j = len(coords)
                coords += geom.coords
                for i in range(j, j + len(geom.coords)-1):
                    segments.append((i, i+1))
            else:
                pass
    elif isinstance(geoms, shapely.geometry.Polygon):
        coords = geoms.exterior.coords
        coords = [c for c in coords]
        for i in range(0, len(coords) - 1):
            segments.append((i, i + 1))

    if output_file_prefix:
        if output_dir is None:
            output_dir = './'

        coords_ofile = os.path.join(output_dir, f'{output_file_prefix}_nodes.json')
        with open(coords_ofile, 'w') as fptr:
            features = []
            id = 1
            for c in coords:
                features.append({
                    "type": "Feature",
                    "id": id,
                    "geometry": shapely.geometry.Point(c[0], c[1]).__geo_interface__,
                    "properties": {}
                })
                id += 1

            fc = {
                "type": "FeatureCollection",
                "features": features
            }
            fptr.write(json.dumps(fc))

        segments_ofile = os.path.join(output_dir, f'{output_file_prefix}_segments.json')
        with open(segments_ofile, 'w') as fptr:
            features = []
            id = 1
            for s in segments:
                p1 = shapely.geometry.Point(coords[s[0]])
                p2 = shapely.geometry.Point(coords[s[1]])

                features.append({
                    "type": "Feature",
                    "id": id,
                    "geometry": shapely.geometry.LineString([p1, p2]).__geo_interface__,
                    "properties": {}
                })
                id += 1

            fc = {
                "type": "FeatureCollection",
                "features": features
            }
            fptr.write(json.dumps(fc))

    return segments, coords


def triangle2shp(tif_file, base_name):
    grid = Grid.from_raster(tif_file)
    crs = grid.crs.srs

    odir = os.path.dirname(tif_file)
    node_file = os.path.join(odir, f'{base_name}.node')
    edge_file = os.path.join(odir, f'{base_name}.edge')
    poly_file = os.path.join(odir, f'{base_name}.poly')
    ele_file = os.path.join(odir, f'{base_name}.ele')

    node_fptr = open(node_file, 'r')
    line = node_fptr.readline()
    num_nodes, num_dims, num_attrs, num_bndry = tuple([int(str(s).strip()) for s in line.split()])
    line = node_fptr.readline()

    nodes = {}
    while line:
        if '#' not in line:
            _idx, _x, _y, _z, _bnd_code = tuple([str(s).strip() for s in line.split()])
            _idx = int(_idx)
            _x = float(_x)
            _y = float(_y)
            _z = float(_z)
            _bnd_code = int(_bnd_code)
            nodes[_idx] = {
                "id": _idx,
                "x": _x,
                "y": _y,
                "bnd": _bnd_code
            }

        line = node_fptr.readline()

    geoms = []
    bnds = []
    for k, v in nodes.items():
        pnt = shapely.geometry.Point(v["x"], v["y"])
        geoms.append(pnt)
        bnds.append(v['bnd'])

    df = pd.DataFrame({'geoms': geoms, 'boundary': bnds})
    gdf = gpd.GeoDataFrame(df, geometry='geoms', crs=crs)
    gdf.to_file(os.path.join(odir, 'points.shp'))

    if os.path.exists(ele_file):
        ele_fptr = open(ele_file, 'r')
        line = ele_fptr.readline()
        num_tris, num_sides, num_bndry = tuple([int(str(s).strip()) for s in line.split()])
        line = ele_fptr.readline()
        eles = {}
        while line:
            if '#' not in line:
                _idx, _node_a, _node_b, _node_c = tuple([int(str(s).strip()) for s in line.split()])
                eles[_idx] = {
                    "id": _idx,
                    "node_a": _node_a,
                    "node_b": _node_b,
                    "node_c": _node_c
                }
            line = ele_fptr.readline()

        geoms = []
        for k, v in eles.items():
            a = v['node_a']
            b = v['node_b']
            c = v['node_c']
            node_a = shapely.geometry.Point(nodes[a]["x"], nodes[a]["y"])
            node_b = shapely.geometry.Point(nodes[b]["x"], nodes[b]["y"])
            node_c = shapely.geometry.Point(nodes[c]["x"], nodes[c]["y"])
            poly = shapely.geometry.Polygon([[p.x, p.y] for p in [node_a, node_c, node_b]])
            geoms.append(poly)

        gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(geoms), crs=crs)
        gdf.to_file(os.path.join(odir, 'tin.shp'))


if __name__ == '__main__':
    wshed_tif = r'/Users/ujjwal/projects/tribs_tin_maker/tin_maker/sample/17/w_17.tif'
    base_name = 'out.1'
    triangle2shp(wshed_tif, base_name)
