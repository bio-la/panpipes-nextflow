#!/usr/bin/env python
import argparse
import os
import sys
import logging
import scanpy as sc
import spatialdata as sd


from funcs.scmethods import run_neighbors_method_choice
from funcs.io import read_yaml
from funcs.scmethods import lsi

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument('--infile',default="spatialdata-preprocessed.zarr",
                    help="file name, format: .zarr")
parser.add_argument('--outfile',default="patialdata-preprocessed.zarr",
                    help="file name, format: .zarr")
parser.add_argument('--neighbor_dict', default=None,
                    help="helps find the correct dimension reduction for sc.pp.neighbors(), default='X_pca'")
parser.add_argument('--n_threads', default=1,help='number of threads available')



args, opt = parser.parse_known_args()
L.info("Running with params: %s", args)

# load the filtering dictionary 
neighbor_dict = args.neighbor_dict
if isinstance(args.neighbor_dict, dict):
    neighbor_dict = args.neighbor_dict
else:
    neighbor_dict = read_yaml(args.neighbor_dict) 

sc.settings.n_jobs = int(args.n_threads)


# read data
L.info("Reading in SpatialData from '%s'" % args.infile)
sdata = sd.read_zarr(args.infile)


 
if neighbor_dict['use_existing'] == "True":
    L.info('Using existing neighbors graph for %s' % args.infile)
    pass
else:
    L.info("Computing new neighbors on %s" % (neighbor_dict['dim_red']))
    if (neighbor_dict['dim_red'] == "X_pca") and ("X_pca" not in sdata["table"].obsm.keys()):
        L.info("X_pca not found, computing it using default parameters")
        sc.tl.pca(sdata["table"])
    opts = dict(method=neighbor_dict['method'],
                n_neighbors=int(neighbor_dict['k']),
                n_pcs=int(neighbor_dict['n_dim_red']),
                metric=neighbor_dict['metric'],
                nthreads=args.n_threads,
                use_rep=neighbor_dict['dim_red'])
    # run command
    run_neighbors_method_choice(sdata["table"],**opts)



L.info("Saving updated SpatialData to '%s'" % args.outfile)
sdata.write(args.outfile)


L.info("Done")
