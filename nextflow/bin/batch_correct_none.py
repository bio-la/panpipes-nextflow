# check for threads number
import multiprocessing 
threads_available = multiprocessing.cpu_count()

# import numpy as np
import pandas as pd
import scanpy as sc
import os
import argparse
from muon import MuData
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice
import muon as mu
import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--input_anndata',
                    default='adata_scaled.h5ad',
                    help='')
parser.add_argument('--modality',
                    default='rna',
                    help='')
parser.add_argument('--dimred',
                    default='PCA',
                    help='which dimred to expect, relevant for ATAC')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_none.csv',
                    help='')
parser.add_argument('--n_threads', default=1, type=int,
                    help="num threads to use for neighbor computations")
parser.add_argument('--integration_col', default='batch')
parser.add_argument('--neighbors_n_pcs',
                    help="n_pcs")
parser.add_argument('--neighbors_method',
                    help="neighbours method, scanpy or hnsw")
parser.add_argument('--neighbors_k',
                    help="neighbors k")
parser.add_argument('--neighbors_metric',
                    help="neighbor metric, e.g. euclidean or cosine")
parser.add_argument('--output_anndata',default=None,
                    help='Path to write the corrected AnnData .h5ad (default: tmp/harmony_scaled_adata_<modality>.h5ad)')


args, opt = parser.parse_known_args()
L.info("Running with params: %s", args)


adata_path = args.input_anndata +"/" + args.modality
if os.path.exists(args.input_anndata):
    L.info("Reading in data from '%s'" % adata_path)
    mdata = mu.read(args.input_anndata)
    mods = list(mdata.mod.keys())
    if args.modality not in mdata.mod:
        L.error("Requested modality '%s' not found. Available modalities: %s",
                args.modality, ", ".join(mods))
        sys.exit(1)
    adata = mdata.mod[args.modality]
    #adata = mu.read(args.input_anndata +"/" + args.modality)
else:
    L.info("missing input anndata")




columns = [x.strip() for x in args.integration_col.split(",")]
if len(columns)>1: 
    comb_columns = "|".join(columns)
    adata.obs[comb_columns] = adata.obs[columns].apply(lambda x: '|'.join(x), axis=1)
    columns += [comb_columns]

# write out batch
# adata.obs[columns].to_csv(os.path.join(os.path.dirname(args.output_csv), 'batch_'+ args.modality +'_mtd.csv'))
if args.dimred == "PCA":
    dimred = "X_pca"
elif args.dimred == "LSI":
    dimred = "X_lsi"


# run neighbours and umap without batch correction
pc_kwargs = {}
if int(args.neighbors_n_pcs) > 0:
    pc_kwargs['use_rep'] = dimred
    pc_kwargs['n_pcs'] = int(args.neighbors_n_pcs)
else:
    # need to push scanpy to use .X if use_rep is None.
    pc_kwargs['use_rep'] = None
    pc_kwargs['n_pcs'] = int(args.neighbors_n_pcs)


if dimred not in adata.obsm:
    L.warning("Dimred '%s' could not be found in adata.obsm. Computing PCA with default parameters." % dimred)
    n_pcs = 50
    if adata.var.shape[0] < n_pcs:
        L.info("You have less features than number of PCs you intend to calculate")
        n_pcs = adata.var.shape[0] - 1
        L.info("Setting n PCS to %i" % int(n_pcs))    
    L.info("Scaling data")
    sc.pp.scale(adata)
    L.info("Computing PCA")
    sc.tl.pca(adata, n_comps=n_pcs, 
                    svd_solver='arpack', 
                    random_state=0) 
    pc_kwargs['use_rep'] = "X_pca"
    pc_kwargs['n_pcs'] = n_pcs

# avoid asking for more pcs than available
if pc_kwargs['use_rep'] is not None:
    avail = adata.obsm[pc_kwargs['use_rep']].shape[1]
    if pc_kwargs['n_pcs'] > avail:
        L.warning(
            "neighbors_n_pcs (%d) > available components in %s (%d); clipping.",
            pc_kwargs['n_pcs'], pc_kwargs['use_rep'], avail
        )
        pc_kwargs['n_pcs'] = avail

L.info("Computing neighbors")
run_neighbors_method_choice(adata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    metric=args.neighbors_metric, 
    nthreads=max([threads_available, 6]), **pc_kwargs)

# Avoid forcing too many threads
# We already define number of cores in nextflow config
# consider adding this change
# L.info("Computing neighbors")
# run_neighbors_method_choice(adata, 
#     method=args.neighbors_method, 
#     n_neighbors=int(args.neighbors_k), 
#     metric=args.neighbors_metric, 
#     nthreads = min(int(args.n_threads), threads_available), **pc_kwargs)

L.info("Computing UMAP")
# For reproducibility
sc.tl.umap(adata, random_state=0)

#write out
L.info("Saving UMAP coordinates to csv file '%s'" % args.output_csv)
umap = pd.DataFrame(adata.obsm['X_umap'], adata.obs.index)
umap.to_csv(args.output_csv)

if args.output_anndata is not None:
    outfile = args.output_anndata
    if os.path.isdir(outfile) or not outfile.endswith('.h5ad'):
        outfile = os.path.join(outfile, f"no_correction_scaled_adata_{args.modality}.h5ad")
else:
    base = os.path.splitext(os.path.basename(args.output_csv))[0]
    outdir = os.path.dirname(args.output_csv)
    outfile = os.path.join(outdir, f"{base}_{args.modality}.h5ad")

#outfiletmp = ("tmp/no_correction_scaled_adata_" + args.modality + ".h5ad" )

L.info("Saving AnnData to '%s'", outfile)
write_anndata(adata, outfile, use_muon=False, modality=args.modality)

L.info("Done")

