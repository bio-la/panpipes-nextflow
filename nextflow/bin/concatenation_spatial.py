#!/usr/bin/env python
'''
Concatenate spatial transcriptomics data
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


import scanpy as sc
import spatialdata as sd

import os
import argparse
import sys
import logging
import glob

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")



sc.settings.verbosity = 3

parser = argparse.ArgumentParser()

parser.add_argument("--input_dirs",
                    default="./tmp/",
                    help="")
parser.add_argument("--output_dir",
                    default="./concatenated.data/",
                    help="")




args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

input_dirs = args.input_dirs.split(" ")

L.info("Reading in all SpatialDatas from %s" % args.input_dirs)
sdatas = []
for file in input_dirs:
    sdatas.append(sd.read_zarr(file))

sdata = sd.concatenate(sdatas, concatenate_tables=True)


L.info("Saving concatenated SpatialData to '%s'" % args.output_dir)
sdata.write(args.output_dir + "concatenated.zarr")

L.info("Done")

