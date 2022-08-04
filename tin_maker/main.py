from process_dem import DemProcessor
from generate_tin import TinGenerator
from utils import get_output_dir, Watershed, run_triangle, triangle2shp
import pprint
import json
import os
import logging


logger = logging.getLogger('tribs_tin_maker')
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

project_name = 'sample'
output_dir = get_output_dir(project_name)

# DEM Processor inputs
padding_percent = 2
flow_accumulation_threshold = 100

# TIN generation inputs
decimation_factor = 0.25
watershed_boundary_buffer_dist = 0.0
stream_boundary_buffer_dist = 0.0


def process_dem(elevation_tif):
    logger.info(f'processing dem {elevation_tif}')
    dem_processor = DemProcessor(dem_path=elevation_tif, tmp_dir=output_dir)
    watersheds = dem_processor.generate_watersheds(padding_percent, flow_accumulation_threshold)
    watersheds_ofile = os.path.join(output_dir, 'watersheds.json')
    with open(watersheds_ofile, 'w') as fptr:
        json_data = []
        for watershed in watersheds:
            json_data.append(watershed.toJson())
        fptr.write(json.dumps(json_data))
    return watersheds_ofile


def generate_tin(watershed_tif, watershed_shp, streams_shp, outlet_csv, odir):
    logger.info(f'processing tin {watershed_tif}')
    tin_generator = TinGenerator(watershed_tif, watershed_shp, streams_shp, outlet_csv)
    nodes_file, edges_file, poly_file = tin_generator.\
        generate_tin(odir, decimation_factor, watershed_boundary_buffer_dist, stream_boundary_buffer_dist)
    run_triangle(poly_file)


if __name__ == '__main__':
    elevation = '/Users/ujjwal/projects/valeriy/clip_dem_utm_12n.tif'
    watersheds_json = process_dem(elevation)
    # watersheds_json = r'/Users/ujjwal/projects/tribs_tin_maker/tin_maker/sample/watersheds.json'
    with open(watersheds_json, 'r') as fptr:
        data = fptr.read()
        watershed_dicts = json.loads(data)

    for str_watershed_dict in watershed_dicts:
        watershed = Watershed()
        watershed_dict = json.loads(str_watershed_dict)
        watershed.from_dict(watershed_dict)
        wshed_number = watershed.id
        odir = watershed.folder
        wshed_tif = watershed.watershed_tif
        catchment_shp = watershed.watershed_shp
        stream_shp = watershed.channel_shp
        outlet_csv = watershed.outlet_csv

        # generate tin
        generate_tin(wshed_tif, catchment_shp, stream_shp, outlet_csv, odir)

        # create shp files for tin and points
        base_name = 'out.1' #triangle runs only once, so output basename will be out.1
        triangle2shp(wshed_tif, base_name)
