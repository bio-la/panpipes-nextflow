import multiprocessing 

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import scvi
import argparse
import os, sys, gc, json, pathlib
import muon as mu
from cgatcore import pipeline as P
import anndata as ad
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
                    help='')
parser.add_argument('--output_csv', default='batch_correction/umap_bc_MultiVI.csv',
                    help='')
parser.add_argument('--integration_col_categorical', default='batch',
                    help='')
parser.add_argument('--integration_col_continuous', default=None,
                    help='')
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
parser.add_argument('--scvi_seed',default=None,
                    help="set explicitly seed to make runs reproducible")

# Replacing YAML params
parser.add_argument('--lowmem', default='true',
                    help="If true, restrict ATAC to ~25k HVFs (like your YAML lowmem)")
# JSON blobs or files (model / training / plan)
parser.add_argument('--model_args_json', default=None)
parser.add_argument('--model_args_json_file', default=None)
parser.add_argument('--training_args_json', default=None)
parser.add_argument('--training_args_json_file', default=None)
parser.add_argument('--training_plan_json', default=None)
parser.add_argument('--training_plan_json_file', default=None)

parser.add_argument('--output_mudata', default=None,
                    help="Path to write .h5mu (dir or filename). Default: tmp/multivi_scaled_adata.h5mu")

# Test params
parser.add_argument('--test_run', default='false',
                    help="If true, override training to a quick test run")
parser.add_argument('--test_max_epochs', type=int, default=10,
                    help="Max epochs to use during --test_run")

args, opt = parser.parse_known_args()
L.info("Running with params: %s", args)

#Helpers
def to_bool(s): return str(s).lower() in ('1','true','t','yes','y')

def _load_json_arg(text=None, file=None):
    if file:
        with open(file) as fh: return json.load(fh)
    if text: return json.loads(text)
    return {}

def _drop_nones(d): return {k: v for k, v in d.items() if v is not None}


# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir

if args.scvi_seed is not None:
    scvi.settings.seed = int(args.scvi_seed)
else:
    scvi.settings.seed = 1492
# load parameters

threads_available = multiprocessing.cpu_count()
#params = pp.io.read_yaml("pipeline.yml")

# ------------------------------------------------------------------
L.info("Reading in MuData from '%s'" % args.scaled_anndata)
mdata = mu.read(args.scaled_anndata)
rna = mdata['rna'].copy()
atac = mdata['atac'].copy()

del mdata

# Optional lowmem ATAC HVF subsetting
if to_bool(args.lowmem):

    if 'highly_variable' in atac.var.columns and atac.var['highly_variable'].any():
        L.info("Subsetting ATAC to existing HVFs")
        atac = atac[:, atac.var.highly_variable].copy()
    L.info("Calculating and subsetting ATAC to top 25k HVFs")
    sc.pp.highly_variable_genes(atac, n_top_genes=25000)
    atac = atac[:, atac.var.highly_variable].copy()

gc.collect()


if rna.shape[0] == atac.shape[0]:
    n=int(rna.shape[0])
else:
    sys.exit("RNA and ATAC have different number of cells, \
        Can't deal with this in this version of MultiVI integration")

gc.collect()

n_genes = len(rna.var_names)
n_regions = len(atac.var_names)


if "raw_counts" in rna.layers.keys():
    L.info("Found raw RNA counts in .layers['raw_counts']")
elif X_is_raw(rna):
    # this means the X layer is already raw and we can make the layer we need
    L.info("Found raw RNA counts in .X. Saving raw RNA counts to .layers['raw_counts']")
    rna.layers["raw_counts"] = rna.X.copy()
else:
    L.error("Could not find raw counts for RNA in .X and .layers['raw_counts']")
    sys.exit("Could not find raw counts for RNA in .X and .layers['raw_counts']")


if "raw_counts" in atac.layers.keys():
    L.info("Found raw ATAC counts in .layers['raw_counts']")
elif X_is_raw(atac):
    # this means the X layer is already raw and we can make the layer we need
    L.info("Found raw ATAC counts in .X. Saving raw ATAC counts to .layers['raw_counts']")
    atac.layers["raw_counts"] = atac.X.copy()
else:
    L.error("Could not find raw counts for ATAC in .X and .layers['raw_counts']")
    sys.exit("Could not find raw counts for ATAC in .X and .layers['raw_counts']")


L.info("Concatenating modalities to comply with multiVI")
# adata_paired = ad.concat([rna, atac], join="outer")
# adata_paired.var = pd.concat([rna.var,atac.var])
if rna.is_view:
    L.info("RNA is view")
    rna = rna.copy()
if atac.is_view:
    L.info("ATAC is view")
    atac = atac.copy()
adata_paired = ad.concat([rna.T, atac.T]).T

rna_cols=rna.obs.columns
atac_cols=atac.obs.columns

rnaobs = rna.obs.copy()
rnaobs.columns= ["rna:"+ x for x in rna_cols]
atacobs = atac.obs.copy()
atacobs.columns= ["atac:"+ x for x in atac_cols]
adata_paired.obs = pd.merge(rnaobs, atacobs, left_index=True, right_index=True)

if "modality" not in adata_paired.obs.columns:
    adata_paired.obs["modality"] = "paired" 



del [rna , atac ]
gc.collect()


L.info("Organizing multiome AnnDatas")
adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired)


# MultiVI integrates by modality, to use batch correction you need a batch covariate to specify in
# categorical_covariate_keys
columns = []
if args.integration_col_categorical is not None :
    #cols = [x.strip() for x in args.integration_col_categorical.split(",")]
    #for cc in cols:
    for cc in [c.strip() for c in args.integration_col_categorical.split(",") if c.strip()]:
        # If already prefixed, keep it — but validate it exists in the merged obs
        if cc.startswith("rna:") or cc.startswith("atac:"):
            if cc in adata_mvi.obs.columns:
                columns.append(cc)
            else:
                raise ValueError(f"Column '{cc}' not found in adata_mvi.obs")
        
        # If not prefixed, detect modality and prefix
        elif cc in rna_cols:
            columns.append(f"rna:{cc}")
        elif cc in atac_cols:
            columns.append(f"atac:{cc}")
        else:
            raise ValueError(f"Column '{cc}' not found in RNA or ATAC obs")

if args.integration_col_continuous is not None :
    if args.integration_col_continuous in rna_cols:
        args.integration_col_continuous = "rna:"+ args.integration_col_continuous
    elif args.integration_col_continuous in atac_cols:
        args.integration_col_continuous = "atac:"+ args.integration_col_continuous


kwargs = {}

if columns is not None:
    print(columns)
    if len(columns) > 1:
        L.info("Using 2 columns to integrate on more variables")
        # bc_batch = "_".join(columns)
        adata_mvi.obs["bc_batch"] = adata_mvi.obs[columns].apply(lambda x: '|'.join(x), axis=1) 
        # make sure that batch is a categorical
        adata_mvi.obs["bc_batch"] = adata_mvi.obs["bc_batch"].astype("category")
    else:
        adata_mvi.obs['bc_batch'] = adata_mvi.obs[columns]
        adata_mvi.obs["bc_batch"] = adata_mvi.obs["bc_batch"].astype("category")

    batch_categories = list(adata_mvi.obs['bc_batch'].unique())
    kwargs['categorical_covariate_keys'] = ["bc_batch"]

if args.integration_col_continuous is not None :
    print(args.integration_col_continuous)
    adata_mvi.obs['bc_batch_continuous'] = adata_mvi.obs[args.integration_col_continuous]
    kwargs['continuous_covariate_keys'] = ["bc_batch_continuous"]


# 1 setup anndata

L.info("Setting up AnnData")
scvi.model.MULTIVI.setup_anndata(
    adata_mvi, 
    batch_key='modality',
    layer =  "raw_counts",
    **kwargs) 
# 2 setup model
model_args     = _drop_nones(_load_json_arg(args.model_args_json,     args.model_args_json_file))
training_args  = _drop_nones(_load_json_arg(args.training_args_json,  args.training_args_json_file))
training_plan  = _drop_nones(_load_json_arg(args.training_plan_json,  args.training_plan_json_file))

# Optional quick test override
if to_bool(args.test_run):
    L.info("Test run enabled: setting max_epochs=%d", int(args.test_max_epochs))
    training_args['max_epochs'] = int(args.test_max_epochs)


L.info("Defining model")
mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=n_genes,
    n_regions=n_regions,
    **model_args
)

#3.train

mvi.view_anndata_setup()

L.info("training args")
print(training_args)


L.info("training plan")
print(training_plan)


L.info("Running multiVI")
mvi.train(**training_args, **training_plan)

mvi.save(os.path.join("batch_correction", "multivi_model"), 
                    anndata=False)

L.info("Plotting ELBO")
plt.plot(mvi.history["elbo_train"], label="train")
plt.plot(mvi.history["elbo_validation"], label="validation")
plt.title("Negative ELBO over training epochs")

plt.legend()
plt.savefig(os.path.join(args.figdir, "multivi_elbo_plot.png"))

L.info("""We support the use of mudata as a general framework for multimodal data
        For this reason, the object we save is not the classical anndata
        you would find in scvitools tutorial.
        We chose to save the learned SC representation in the 
        Mudata obsm slot, and any single modality processing in
        its own modality slot
            """)

mdata = mu.read(args.scaled_anndata)
L.info("Extracting latent space and saving latent to X_MultiVI")
mdata.obsm["X_MultiVI"] = mvi.get_latent_representation()


if int(args.neighbors_n_pcs) > mdata.obsm['X_MultiVI'].shape[1]:
    L.warn(f"N PCs is larger than X_MultiVI dimensions, reducing n PCs to  {mdata.obsm['X_MultiVI'].shape[1] -1}")
n_pcs= min(int(args.neighbors_n_pcs), mdata.obsm['X_MultiVI'].shape[1]-1)

L.info("Computing neighbors")
run_neighbors_method_choice(mdata, 
    method=args.neighbors_method, 
    n_neighbors=int(args.neighbors_k), 
    n_pcs=n_pcs, 
    metric=args.neighbors_metric, 
    use_rep='X_MultiVI',
    nthreads=max([threads_available, 6]))
L.info("Computing UMAP")
sc.tl.umap(mdata, min_dist=0.4)
L.info("Computing Leiden clustering")
sc.tl.leiden(mdata, key_added="leiden_multiVI")

L.info("Saving UMAP coordinates to csv file '%s" % args.output_csv)
umap = pd.DataFrame(mdata.obsm['X_umap'], mdata.obs.index)
umap.to_csv(args.output_csv)


if args.output_mudata:
    out_path = args.output_mudata
    if os.path.isdir(out_path) or not out_path.endswith('.h5mu'):
        out_path = os.path.join(out_path, "multivi_scaled_adata.h5mu")
else:
    out_path = "tmp/multivi_scaled_adata.h5mu"

out_dir = os.path.dirname(out_path)
if out_dir:
    os.makedirs(out_dir, exist_ok=True)

L.info("Saving MuData to '%s'", out_path)
mdata.write(out_path)


L.info("Done")

