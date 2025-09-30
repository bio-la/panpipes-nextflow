import multiprocessing 

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import gc
import json
import muon as mu
from cgatcore import pipeline as P
import panpipes.funcs as pp
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice, X_is_raw

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# load arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--scaled_anndata',
                    default='adata_scaled.h5ad',
                    help='a preprocessed mudata/anndata')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_MultiVI.csv',
                    help='')
parser.add_argument('--use_gpu', default=False,
                    help='')
parser.add_argument('--integration_col_categorical', default='batch',
                    help='')
parser.add_argument('--n_factors', default="",
                    help='number of factors to train the model with - optional')
parser.add_argument('--n_iterations', default='',
                    help='upper limit on the number of iterations')
parser.add_argument('--convergence_mode', default='fast',
                    help='fast, medium or slow')
parser.add_argument('--save_parameters', default=None,
                    help='whether to save parameters models')
parser.add_argument('--outfile_model', default="",
                    help="if args.save_parameters true, path to HDF5 file to store the model")
parser.add_argument('--figdir', default='./figures',
                    help='')
parser.add_argument('--neighbors_n_pcs',
                    help="n_pcs", default=50)
parser.add_argument('--neighbors_method', default="scanpy",
                    help="neighbours method, scanpy or hnsw")
parser.add_argument('--neighbors_k', default=30,
                    help="neighbors k")
parser.add_argument('--neighbors_metric',default="euclidean",
                    help="neighbor metric, e.g. euclidean or cosine")

# Replace YAML
parser.add_argument('--mofa_args_json', default=None)
parser.add_argument('--mofa_args_json_file', default=None)
parser.add_argument('--modalities', default=None,
                    help="Comma-separated list of modalities to use (default: all)")
parser.add_argument('--output_mudata', default=None,
                    help="Path to write .h5mu (dir or filename). Default: tmp/mofa_scaled_adata.h5mu")
parser.add_argument('--filter_by_hvg', default='true',
                    help="If true, require/use 'highly_variable' per modality (default: true)")


args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)
# Helpers
def to_bool(s):
    if isinstance(s, bool):
        return s
    if s is None:
        return False
    return str(s).lower() in ('1','true','t','yes','y')

def _load_json_arg(text=None, file=None):
    if file:
        with open(file) as fh:
            return json.load(fh)
    if text:
        return json.loads(text)
    return {}

def _drop_nones(d):
    return {k: v for k, v in d.items() if v is not None}

# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir


# load parameters

threads_available = multiprocessing.cpu_count()
#params = pp.io.read_yaml("pipeline.yml")

L.info("Reading in MuData from '%s'" % args.scaled_anndata)
mdata = mu.read(args.scaled_anndata)

# if params['multimodal']['mofa']['modalities'] is not None:
#     modalities= params['multimodal']['mofa']['modalities']
#     modalities = [x.strip() for x in modalities.split(",")]
#     L.info(f"Using modalities :{modalities}")
#     removed_mods = None
#     if all(x in modalities for x in mdata.mod.keys()):
#         tmp = mdata.copy()
#         L.info('Using all modalities')
#     else:
#         tmp = mdata.copy()
#         removed_mods = list(set(mdata.mod.keys()) - set(modalities))
#         L.info(f"Removing modalities {removed_mods}")
#         for rmod in removed_mods:
#             del tmp.mod[rmod]
# else:
#     L.warning("""No modalities were specified, therefore all available modalities will be used.
#                 This may be problematic if repertoire is present as modality.""")
#     removed_mods = None
#     tmp = mdata.copy()  

if args.modalities:
    modalities = [x.strip() for x in str(args.modalities).split(",") if x.strip()]
    L.info(f"Using modalities: {modalities}")
    removed_mods = None
    if all(x in mdata.mod.keys() for x in modalities):
        tmp = mdata.copy()
        # drop non-selected mods from tmp
        for rmod in list(set(mdata.mod.keys()) - set(modalities)):
            del tmp.mod[rmod]
        if len(set(mdata.mod.keys()) - set(modalities)) > 0:
            removed_mods = list(set(mdata.mod.keys()) - set(modalities))
            L.info(f"Removed modalities from tmp: {removed_mods}")
    else:
        # remove any not requested
        tmp = mdata.copy()
        removed_mods = list(set(mdata.mod.keys()) - set(modalities))
        L.info(f"Removing modalities {removed_mods}")
        for rmod in removed_mods:
            if rmod in tmp.mod:
                del tmp.mod[rmod]
else:
    L.warning("""No modalities were specified, therefore all available modalities will be used.
                This may be problematic if repertoire is present as modality.""")
    removed_mods = None
    tmp = mdata.copy()


L.info("Intersecting modality obs before running mofa")
mu.pp.intersect_obs(tmp)

#mofa_kwargs={}
mofa_kwargs = _drop_nones(_load_json_arg(args.mofa_args_json, args.mofa_args_json_file))
mofa_kwargs.pop('neighbors', None)
mofa_kwargs.pop('modalities', None)


#expected args:
# n_factors: 10
# n_iterations: 1000
# convergence_mode: fast, medium, slow
# save_parameters: False
# #if save_parameters True, set the following, otherwise leave blank
# outfile_model:

# mofa_kwargs = params['multimodal']['mofa']
# del mofa_kwargs['modalities']

# if mofa_kwargs['filter_by_hvg']:
#     mofa_kwargs['use_var'] = "highly_variable"
#     del mofa_kwargs['filter_by_hvg']
#     for mod in tmp.mod.keys():
#         if "highly_variable" not in tmp[mod].var.columns:
#             tmp[mod].var["highly_variable"] = True
        
#     tmp.update()

if 'filter_by_hvg' in mofa_kwargs:
    filter_by_hvg = to_bool(mofa_kwargs.pop('filter_by_hvg'))
else:
    filter_by_hvg = to_bool(args.filter_by_hvg)

if filter_by_hvg:
    for mod in list(tmp.mod.keys()):
        if 'highly_variable' not in tmp[mod].var.columns:
            tmp[mod].var['highly_variable'] = True
        hv_mask = tmp[mod].var['highly_variable'].values
        tmp.mod[mod] = tmp.mod[mod][:, hv_mask].copy()
    tmp.update()
    
mofa_kwargs.pop('use_var', None)

_ALLOWED = {
    'n_factors', 'n_iterations', 'convergence_mode',
    'save_parameters', 'outfile', 'use_var', 'groups_label', 'gpu_mode'
}
mofa_kwargs = {k: v for k, v in mofa_kwargs.items() if k in _ALLOWED}


if args.integration_col_categorical is not None:
    if args.integration_col_categorical in tmp.obs.columns:
        mofa_kwargs['groups_label'] = args.integration_col_categorical

# default is to read yaml and parse directly to kwargs.
# if the defaults expected params are parsed by the script in some other way
# they will overwrite the initial reading of the yml

# if mofa_kwargs['save_parameters'] is None:
#     if args.save_parameters is not None:
#         mofa_kwargs['save_parameters'] = check_for_bool(args.save_parameters)
#         if args.outfile_model is not None:
#             mofa_kwargs['outfile'] = args.outfile_model

def _nonempty(x):
    return x is not None and str(x).strip() != ""

# save_parameters / outfile
sp = mofa_kwargs.get('save_parameters', None)
if sp is None and _nonempty(args.save_parameters):
    mofa_kwargs['save_parameters'] = check_for_bool(args.save_parameters)
    if _nonempty(args.outfile_model):
        mofa_kwargs['outfile'] = args.outfile_model

# added str(args.xxx).strip() because Nextflow/Groovy often passes “unset”
# params as empty strings, not None

# n_factors
if mofa_kwargs.get('n_factors', None) is None and _nonempty(args.n_factors):
    mofa_kwargs['n_factors'] = int(args.n_factors)

# n_iterations
if mofa_kwargs.get('n_iterations', None) is None and _nonempty(args.n_iterations):
    mofa_kwargs['n_iterations'] = int(args.n_iterations)

# convergence_mode
if mofa_kwargs.get('convergence_mode', None) is None and _nonempty(args.convergence_mode):
    mofa_kwargs['convergence_mode'] = str(args.convergence_mode)

# are we using the gpu?
# gpu
if _nonempty(args.use_gpu):
    mofa_kwargs['gpu_mode'] = to_bool(args.use_gpu)


L.info("Running mofa")
mu.tl.mofa(tmp, **mofa_kwargs)
#This adds X_mofa embeddings to the .obsm slot of the MuData object
#write the discovered latent rep to the original mudata object
L.info("Adding X_mofa to .obsm of MuData")
mdata.obsm['X_mofa'] = tmp.obsm['X_mofa'].copy()

if int(args.neighbors_n_pcs) > mdata.obsm['X_mofa'].shape[1]:
    L.warning(f"N PCs is larger than X_mofa dimensions, reducing n PCs to  {mdata.obsm['X_mofa'].shape[1]-1}")
n_pcs= min(int(args.neighbors_n_pcs), mdata.obsm['X_mofa'].shape[1]-1)

L.info("Computing neighbors")
run_neighbors_method_choice(mdata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=n_pcs, #this should be the # rows of var, not obs ???????
    metric=args.neighbors_metric, 
    use_rep='X_mofa',
    nthreads=max([threads_available, 6]))

L.info("Computing UMAP")
sc.tl.umap(mdata, min_dist=0.4)
L.info("Computing Leiden clustering")
sc.tl.leiden(mdata,  key_added="leiden_mofa")

L.info("Saving UMAP coordinates to csv file '%s" % args.output_csv)
umap = pd.DataFrame(mdata.obsm['X_umap'], mdata.obs.index)
umap.to_csv(args.output_csv)

if removed_mods is not None:
    for rmd in removed_mods:
        tmp.mod[rmd] = mdata.mod[rmd].copy()

if args.output_mudata:
    out_path = args.output_mudata
    if os.path.isdir(out_path) or not out_path.endswith('.h5mu'):
        out_path = os.path.join(out_path, "mofa_scaled_adata.h5mu")
else:
    out_path = "tmp/mofa_scaled_adata.h5mu"

out_dir = os.path.dirname(out_path)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)

L.info("Saving MuData to 'tmp/mofa_scaled_adata.h5mu'")
mdata.write(out_path)

L.info("Done")

