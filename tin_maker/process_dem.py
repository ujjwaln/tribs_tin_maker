import numpy as np
import os
import geopandas as gpd
from pysheds.grid import Grid
from pysheds import io
from geojson import LineString
import fiona
import fiona.crs
from tin_maker.utils import geojson2shapely, Watershed


class DemProcessor(object):

    def __init__(self, dem_path, tmp_dir):
        self.dem_path = dem_path
        self.tmp_dir = tmp_dir

    def generate_watersheds(self, padding_percent=2, flow_accumulation_threshold=50):
        """
        generate subwatersheds in rectangular DEM. Works by computing channels based on flow accumulation
        threshold and then intersecting channels with the DEM rectangluar extent boundary to locate outlets.
        For each outlet, the upstream watershed is computed.

        :param elevation: GeoTiff raster digital elevation model
        :param padding_percent: When computing outlets, shrink raster extent by this % value
        :param flow_accumulation_threshold: Flow acc. threshold when creating channels (# of cells)
        :return:
        """

        # initialize the grid
        elevation = self.dem_path
        grid = Grid.from_raster(elevation)

        # check coordinate system is cartesian
        assert (grid.crs.crs.coordinate_system.name == 'cartesian')

        # read the raster
        dem = grid.read_raster(elevation)
        xmin, xmax, ymin, ymax = dem.extent

        padding = (xmax - xmin) * padding_percent / 100.
        xmin = xmin + padding
        xmax = xmax - padding
        ymax = ymax - padding
        ymin = ymin + padding

        line_top = LineString(coordinates=[(xmin, ymax), (xmax, ymax)])
        line_top = geojson2shapely(line_top)

        line_bottom = LineString(coordinates=[(xmin, ymin), (xmax, ymin)])
        line_bottom = geojson2shapely(line_bottom)

        line_left = LineString(coordinates=[(xmin, ymin), (xmin, ymax)])
        line_left = geojson2shapely(line_left)

        line_right = LineString(coordinates=[(xmax, ymin), (xmax, ymax)])
        line_right = geojson2shapely(line_right)

        # Condition DEM
        # ----------------------
        # Fill pits in DEM
        pit_filled_dem = grid.fill_pits(dem)

        # Fill depressions in DEM
        flooded_dem = grid.fill_depressions(pit_filled_dem)

        # Resolve flats in DEM
        inflated_dem = grid.resolve_flats(flooded_dem)

        # save for resampling elevations
        ofile = os.path.join(os.path.dirname(elevation),
                             '%s_fill.tif' % (os.path.splitext(os.path.basename(elevation))[0]))
        grid.to_raster(inflated_dem, ofile, blockxsize=16, blockysize=16)

        # Determine D8 flow directions from DEM
        # ----------------------
        # Specify directional mapping
        dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

        # Compute flow directions
        # -------------------------------------
        fdir = grid.flowdir(inflated_dem, dirmap=dirmap)

        # Calculate flow accumulation
        # --------------------------
        acc = grid.accumulation(fdir, dirmap=dirmap)

        # Extract river network
        # ---------------------
        # branches is type geojson.features.FeatureCollection
        branches = grid.extract_river_network(fdir, acc > flow_accumulation_threshold, dirmap=dirmap)

        # Write shapefile
        schema = {
            'geometry': 'LineString',
            'properties': {'LABEL': 'float:16'}
        }
        channels_file = os.path.join(self.tmp_dir, 'channels.shp')
        with fiona.open(channels_file, 'w',
                        driver='ESRI Shapefile',
                        crs=grid.crs.srs,
                        schema=schema) as c:
            j = 0
            for f in branches.features:
                rec = {'geometry': f.geometry, 'properties': {'LABEL': str(j)}, 'id': str(j)}
                c.write(rec)
                j += 1

        # extract outlets
        channels_gdf = gpd.read_file(channels_file)
        outlets = []
        for branch in branches['features']:
            shapely_branch = geojson2shapely(branch['geometry'])
            for line in [line_bottom, line_right, line_left, line_top]:
                intersection = line.intersection(shapely_branch)
                if intersection:
                    outlets.append(intersection.centroid)

        i = 0
        watersheds = []
        for outlet in outlets:
            try:
                catch = grid.catchment(x=outlet.x, y=outlet.y, fdir=fdir, xytype='coordinate')
                # save
                inflated_dem.mask = catch
                _odir = os.path.join(self.tmp_dir, str(i))
                os.makedirs(_odir, exist_ok=True)

                ofile = os.path.join(_odir, f'w_{i}.tif')
                if os.path.exists(ofile):
                    os.remove(ofile)
                io.to_raster(inflated_dem, ofile, dtype=np.float)

                # Clip to catchment
                grid.clip_to(catch)

                # Create view
                catch_view = grid.view(catch, dtype=np.uint8)

                # Create a vector representation of the catchment mask
                shapes = grid.polygonize(catch_view)

                # Specify schema
                schema = {
                    'geometry': 'Polygon',
                    'properties': {'LABEL': 'float:16'}
                }

                # Write shapefile
                catchment_file = os.path.join(_odir, f'catchment_{i}.shp')
                with fiona.open(catchment_file, 'w',
                                driver='ESRI Shapefile',
                                crs=grid.crs.srs,
                                schema=schema) as c:
                    j = 0
                    for shape, value in shapes:
                        rec = {'geometry': shape, 'properties': {'LABEL': str(value)}, 'id': str(j)}
                        c.write(rec)
                        j += 1

                grid.viewfinder = dem.viewfinder

                # clip channels
                catchment_gdf = gpd.read_file(catchment_file)
                catchment_geom = catchment_gdf.geometry[0]

                cell_size = (grid.affine[0] + -1 * grid.affine[4]) * 0.5

                # first find intersection
                channels_clipped_gdf = channels_gdf.intersection(catchment_geom)

                # now ensure intersected lines are completely within
                mask = channels_clipped_gdf.geometry.intersects(catchment_geom.buffer(-0.1 * cell_size)) == True
                channels_clipped_gdf = channels_clipped_gdf.loc[mask]

                channels_clipped_file = os.path.join(_odir, f'channels_{i}.shp')
                channels_clipped_gdf.to_file(channels_clipped_file, crs=grid.crs.srs)

                # save outlet location
                outlet_file = os.path.join(_odir, f'outlet_{i}.csv')
                with open(outlet_file, 'w') as ofptr:
                    ofptr.write(f'{outlet.x}, {outlet.y}')

                watershed = Watershed(ofile, catchment_file, channels_clipped_file, outlet_file, catchment_geom.area, i)
                watersheds.append(watershed)

            except ValueError as ex:
                # pour point out of bounds
                print(str(ex))
                pass

            i += 1

        return watersheds

