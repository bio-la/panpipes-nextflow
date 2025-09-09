import argparse, os, sys, json, gc, multiprocessing, pathlib

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import muon as mu
from cgatcore import pipeline as P
import panpipes.funcs as pp
from panpipes.funcs.processing import check_for_bool
from panpipes.funcs.io import read_anndata, write_anndata
from panpipes.funcs.scmethods import run_neighbors_method_choice, X_is_raw

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
parser.add_argument('--output_csv', default='batch_correction/umap_bc_wnn.csv',
                    help='')
parser.add_argument('--figdir', default='./figures',
                    help='')
#Gloabal params
parser.add_argument('--n_neighbors',
                    help="number of neighbors", default=50)
parser.add_argument('--n_bandwidth_neighbors', default="",
                    help="")
parser.add_argument('--n_multineighbors', default=30,
                    help="neighbors k")
parser.add_argument('--metric',default="euclidean",
                    help="neighbor metric, e.g. euclidean or cosine")
parser.add_argument('--low_memory',default=True,
                    help="set to True by default if cells in dataset >50k")

# Modalities to use in WNN (must exist in object)
parser.add_argument('--modalities', default='rna,prot,atac', help='Comma-separated modalities to include, e.g. "rna,prot"')

# RNA
parser.add_argument('--rna_neighbors_method', default='scanpy')
parser.add_argument('--rna_neighbors_metric', default='euclidean')
parser.add_argument('--rna_neighbors_n_pcs', type=int, default=30)
parser.add_argument('--rna_neighbors_k', type=int, default=30)

# Prot
parser.add_argument('--prot_neighbors_method', default='scanpy')
parser.add_argument('--prot_neighbors_metric', default='euclidean')
parser.add_argument('--prot_neighbors_n_pcs', type=int, default=30)
parser.add_argument('--prot_neighbors_k', type=int, default=30)

# ATAC
parser.add_argument('--atac_neighbors_method', default='scanpy')
parser.add_argument('--atac_neighbors_metric', default='euclidean')
parser.add_argument('--atac_neighbors_n_pcs', type=int, default=30)
parser.add_argument('--atac_neighbors_k', type=int, default=30)
parser.add_argument('--atac_dimred', default=None, help='Preferred embedding key for ATAC (e.g., "X_lsi")')


# Batch-corrected choices (replace YAML):
# Pass mapping like {"rna": "bbknn", "prot": "harmony", "atac": null}
# This determines which obsm to expect: "X_<METHOD>" (with "X_scVI" capitalization fix)
parser.add_argument('--batch_corrected_json', default=None, help='JSON dict modality->method or null, e.g. {"rna":"bbknn","prot":"harmony","atac":null}')
parser.add_argument('--batch_corrected_json_file', default=None)

# Optional map of precomputed per-modality anndata paths (if embeddings are not present in current MuData):
# e.g. {"rna": "tmp/harmony_scaled_adata_rna.h5ad", "prot": "tmp/bbknn_scaled_adata_prot.h5ad"}
parser.add_argument('--precomputed_anndata_json', default=None)
parser.add_argument('--precomputed_anndata_json_file', default=None)

parser.add_argument('--output_mudata', default=None, help='Output MuData path (.h5mu). If dir or wrong suffix, saves to <dir>/wnn_scaled_adata.h5mu')

args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)
# scanpy settings
sc.set_figure_params(facecolor="white")
sc.settings.autoshow = False
sc.settings.figdir = args.figdir

#params = pp.io.read_yaml("pipeline.yml")

# Helpers
def to_bool(x):
    return str(x).lower() in ("1", "true", "t", "yes", "y")

def _load_json_arg(text=None, file=None):
    if file:
        with open(file) as fh:
            return json.load(fh)
    if text:
        return json.loads(text)
    return {}

def _drop_nones(d):
    return {k: v for k, v in (d or {}).items() if v is not None}

def _ensure_dir(p):
    if p and not os.path.exists(p):
        os.makedirs(p, exist_ok=True)

def _first_present(*xs):
    for x in xs:
        if x is not None:
            return x
    return None

# ensure requested representation exists (PCA/LSI)
def _ensure_representation(adata, mod, repuse, npcs):
    """
    Ensure 'repuse' exists in adata.obsm.
    If 'X_lsi' requested but missing, fall back to 'X_pca' and compute it.
    Returns the final repuse to use.
    """
    if repuse == "X_pca":
        if "X_pca" not in adata.obsm:
            L.info("[%s] Computing PCA because X_pca is missing", mod)
            sc.tl.pca(adata, n_comps=max(int(npcs), 10))
    elif repuse == "X_lsi":
        if "X_lsi" not in adata.obsm:
            L.warning("[%s] Requested X_lsi but not present; falling back to PCA", mod)
            if "X_pca" not in adata.obsm:
                sc.tl.pca(adata, n_comps=max(int(npcs), 10))
            repuse = "X_pca"
    return repuse

threads_available = multiprocessing.cpu_count()


# if params['multimodal']['WNN']['modalities'] is not None:
#     modalities= params['multimodal']['WNN']['modalities']
#     modalities = [x.strip() for x in modalities.split(",")]
#     L.info(f"Using modalities :{modalities}")

# wnn_params_bc = params['multimodal']['WNN']['batch_corrected'] 
# if modalities is not None:
#     wnn_params_bc= {k: wnn_params_bc[k] for k in wnn_params_bc.keys() & modalities}
# L.info( wnn_params_bc )

L.info("Reading in MuData from '%s'" % args.scaled_anndata)
mdata = mu.read(args.scaled_anndata)

# Modalities to use from args
mods_requested = [m.strip() for m in str(args.modalities).split(",") if m.strip()]
mods_present   = list(mdata.mod.keys())
mods           = [m for m in mods_requested if m in mods_present]
if not mods:
    raise ValueError(...)
L.info("Using modalities for WNN: %s", mods)

# if all(x in modalities for x in mdata.mod.keys()):
#     tmp = mdata.copy()
#     removed_mods = None
# else:
#     tmp = mdata.copy()
#     removed_mods = list(set(mdata.mod.keys()) - set(modalities))
#     L.info(f"Removing modalities {removed_mods}")
#     for rmod in removed_mods:
#         del tmp.mod[rmod]

# Copy and drop unused modalities
tmp = mdata.copy()
removed_mods = list(set(mods_present) - set(mods))
if removed_mods:
    L.info("Removing modalities (not used in this run): %s", removed_mods)
    for r in removed_mods:
        del tmp.mod[r]


L.info("Intersecting modality obs before running WNN")
mu.pp.intersect_obs(tmp)

# batch-corrected mapping
bc_map = _drop_nones(_load_json_arg(args.batch_corrected_json,
                                    args.batch_corrected_json_file))
bc_map = {k: bc_map.get(k, None) for k in mods} if bc_map else {k: None for k in mods}
L.info("Batch-corrected choices per modality: %s", bc_map)

precomp_map = _drop_nones(_load_json_arg(args.precomputed_anndata_json, args.precomputed_anndata_json_file)) or {}


#one could also check of obsp is not empty and use precomputed connectivities but i assume that if batch corrected the obsp is populated, if not, calc on the flight on pca
# dict_graph = {}
# for x in wnn_params_bc.keys():
#     dict_graph[x] = {}
#     if wnn_params_bc[x] is not None:
#         dict_graph[x]["obsm"] = "X_" + wnn_params_bc[x]
#         dict_graph[x]["anndata"] = "tmp/" + wnn_params_bc[x].lower() + "_scaled_adata_" + x +".h5ad"
#     else: 
#         dict_graph[x]["obsm"] = None

# L.info(dict_graph)

# if dict_graph["rna"]["obsm"] == "X_scvi":
#     dict_graph["rna"]["obsm"] = "X_scVI"

dict_graph = {}
for mod in mods:
    method = bc_map.get(mod, None)
    dict_graph[mod] = {"obsm": None, "anndata": None}
    if method is not None:
        # scVI capitalisation fix
        obsm_key = "X_scVI" if method.lower() == "scvi" else f"X_{method}"
        dict_graph[mod]["obsm"] = obsm_key
        # prefer explicit path if given, else use your old tmp/<method>_scaled_adata_<mod>.h5ad default
        dict_graph[mod]["anndata"] = f"tmp/{method.lower()}_scaled_adata_{mod}.h5ad"

L.info("dict_graph (per-modality embedding/backup path): %s", dict_graph)

# ... earlier: mdata read, mods selection, tmp = mdata.copy(), intersect_obs(tmp),
# bc_map -> dict_graph, etc. (unchanged ordering as your original)

threads_available = multiprocessing.cpu_count()

for kmod in dict_graph.keys():
    L.info(kmod)
    # replace params['multimodal']['WNN']['knn'][kmod] with CLI flags
    if kmod == "rna":
        pkmod = {
            'method': args.rna_neighbors_method,
            'metric': args.rna_neighbors_metric,
            'npcs'  : int(args.rna_neighbors_n_pcs),
            'k'     : int(args.rna_neighbors_k),
        }
    elif kmod == "prot":
        pkmod = {
            'method': args.prot_neighbors_method,
            'metric': args.prot_neighbors_metric,
            'npcs'  : int(args.prot_neighbors_n_pcs),
            'k'     : int(args.prot_neighbors_k),
        }
    elif kmod == "atac":
        pkmod = {
            'method': args.atac_neighbors_method,
            'metric': args.atac_neighbors_metric,
            'npcs'  : int(args.atac_neighbors_n_pcs),
            'k'     : int(args.atac_neighbors_k),
        }
    else:
        raise ValueError(f"Unknown modality {kmod}")

    repuse = None
    if dict_graph[kmod]["obsm"] is not None:
        if dict_graph[kmod]["obsm"] not in tmp.mod[kmod].obsm.keys():
            L.info("Provided mdata doesn't have the desired obsm, just checking if it's bbknn you want.")
            if dict_graph[kmod]["obsm"] == "X_bbknn":
                if len(tmp.mod[kmod].obsp.keys()) > 0 and "neighbors" in tmp.mod[kmod].uns.keys():
                    L.info("Populated obsp slot found. Assuming it's bbknn")
                else:
                    if dict_graph[kmod]["anndata"] is not None and os.path.exists(dict_graph[kmod]["anndata"]):
                        L.info("Reading precomputed connectivities for bbknn")
                        adata = mu.read(dict_graph[kmod]["anndata"])
                        tmp.mod[kmod].obsp = adata.obsp.copy()
                        if dict_graph[kmod]["obsm"] in adata.obsm:
                            tmp.mod[kmod].obsm[dict_graph[kmod]["obsm"]] = adata.obsm[dict_graph[kmod]["obsm"]].copy()
                        if "neighbors" in adata.uns:
                            tmp.mod[kmod].uns["neighbors"]= adata.uns["neighbors"].copy()
                        tmp.update()
            else:
                if dict_graph[kmod]["anndata"] is not None and os.path.exists(dict_graph[kmod]["anndata"]):
                    L.info("Provided mdata doesn't have the desired obsm. Reading the batch-corrected data from another stored object")
                    adata = mu.read(dict_graph[kmod]["anndata"])
                    L.debug(kmod + " object")
                    L.debug(adata)
                    if dict_graph[kmod]["obsm"] in adata.obsm:
                        tmp.mod[kmod].obsm[dict_graph[kmod]["obsm"]] = adata.obsm[dict_graph[kmod]["obsm"]].copy()
                    tmp.mod[kmod].obsp = adata.obsp.copy()
                    if "neighbors" in adata.uns:
                        tmp.mod[kmod].uns['neighbors'] = adata.uns['neighbors'].copy()
                    tmp.update()
                    repuse = dict_graph[kmod]["obsm"]
                else:
                    L.warning("Could not find the desired obsm and the anndata slot is empty, will calculate on the fly")
                    if kmod == "atac":
                        if "X_lsi" in tmp.mod[kmod].obsm.keys():
                            repuse = "X_lsi"
                        else:
                            repuse = "X_pca"
                        L.info("Falling back on %s" % (repuse))
                    L.info("Computing neighbours")
                    if repuse != "X_bbknn":
                        run_neighbors_method_choice(
                            tmp.mod[kmod],
                            method=pkmod['method'],
                            n_neighbors=int(pkmod['k']),
                            n_pcs=min(int(pkmod['npcs']), tmp.mod[kmod].var.shape[0]-1),
                            metric=pkmod['metric'],
                            use_rep=repuse,
                            nthreads=max([threads_available, 6])
                        )
        else:
            L.info("Using %s" % (dict_graph[kmod]["obsm"]))
    else:
        L.warning("Could not find the desired obsm and the anndata slot is empty, will calculate on the fly")
        repuse = "X_pca"
        if kmod == "atac":
            if "X_lsi" in tmp.mod[kmod].obsm.keys():
                repuse = "X_lsi"
            else:
                repuse = "X_pca"
            L.info("falling back on %s" % (repuse))
        L.info("Computing neighbours")
        repuse = _ensure_representation(tmp.mod[kmod], kmod, repuse, pkmod['npcs'])
        run_neighbors_method_choice(
            tmp.mod[kmod],
            method=pkmod['method'],
            n_neighbors=int(pkmod['k']),
            n_pcs=min(int(pkmod['npcs']), tmp.mod[kmod].var.shape[0]-1),
            metric=pkmod['metric'],
            use_rep=repuse,
            nthreads=max([threads_available, 6])
        )

tmp.update()


# for kmod in dict_graph.keys():
#     L.info(kmod)
#     pkmod=params['multimodal']['WNN']['knn'][kmod]
    
#     if dict_graph[kmod]["obsm"] is not None:
#         if dict_graph[kmod]["obsm"] not in tmp.mod[kmod].obsm.keys():
#             L.info("Provided mdata doesn't have the desired obsm, just checking if it's bbknn you want.")
#             if dict_graph[kmod]["obsm"] == "X_bbknn":
#                 if len(tmp.mod[kmod].obsp.keys()) > 0 and "neighbors" in tmp.mod[kmod].uns.keys() : 
#                      L.info("Populated obsp slot found. Assuming it's bbknn")
#                 else:
#                     if dict_graph[kmod]["anndata"] is not None:
#                         L.info("Reading precomputed connectivities for bbknn")
#                         adata = mu.read(dict_graph[kmod]["anndata"])
#                         tmp.mod[kmod].obsp = adata.obsp.copy()
#                         tmp.mod[kmod].obsm[dict_graph[kmod]["obsm"]] = adata.obsm[dict_graph[kmod]["obsm"]].copy()
#                         tmp.mod[kmod].uns["neighbors"]= adata.uns["neighbors"].copy()
#                         tmp.update()
#             else:
#                 if dict_graph[kmod]["anndata"] is not None:
#                     L.info("Provided mdata doesn't have the desired obsm. reading the batch corrected data from another stored object")
#                     adata = mu.read(dict_graph[kmod]["anndata"])
#                     L.debug(kmod + "object")
#                     L.debug(adata)
#                     tmp.mod[kmod].obsm[dict_graph[kmod]["obsm"]] = adata.obsm[dict_graph[kmod]["obsm"]].copy()
#                     tmp.mod[kmod].obsp = adata.obsp.copy()
#                     tmp.mod[kmod].uns['neighbors'] = adata.uns['neighbors'].copy()
#                     tmp.update()
#                     repuse = dict_graph[kmod]["obsm"] 
#                 else:
#                     L.warning("Could not find the desired obsm and the anndata slot is empty, will calculate on the flight")
#                     if kmod =="atac":
#                         if "X_lsi" in tmp.mod[kmod].obsm.keys():
#                             repuse = "X_lsi"
#                         else:
#                             repuse = "X_pca"
#                         L.info("Falling back on %s" %(repuse) )
            
#                     L.info("Computing neighbours")
#                     if repuse != "X_bbknn":
#                         run_neighbors_method_choice(tmp.mod[kmod], 
#                             method=pkmod['method'], 
#                             n_neighbors=int(pkmod['k']), 
#                             n_pcs=min(int(pkmod['npcs']), mdata.var.shape[0]-1), #this should be the # rows of var, not obs
#                             metric=pkmod['metric'], 
#                             #does this throw an error if no PCA for any single mod is stored?
#                             use_rep=repuse,
#                             nthreads=max([threads_available, 6]))
#         else:
#             L.info("Using %s" %(dict_graph[kmod]["obsm"]))
#     else:
#         L.warning("Could not find the desired obsm and the anndata slot is empty, will calculate on the flight")
#         repuse ="X_pca"
#         if kmod =="atac":
#             if "X_lsi" in tmp.mod[kmod].obsm.keys():
#                 repuse = "X_lsi"
#             else:
#                 repuse = "X_pca"
#             L.info("falling back on %s" %(repuse) )

#         L.info("Computing neighbours")
#         run_neighbors_method_choice(tmp.mod[kmod], 
#             method=pkmod['method'], 
#             n_neighbors=int(pkmod['k']), 
#             n_pcs=min(int(pkmod['npcs']), mdata.var.shape[0]-1), #this should be the # rows of var, not obs
#             metric=pkmod['metric'], 
#             #does this throw an error if no PCA for any single mod is stored?
#             use_rep=repuse,
#             nthreads=max([threads_available, 6]))

# tmp.update()

low_mem = args.low_memory
if isinstance(low_mem, str):
    low_mem = True if low_mem.lower() in ("1","true","t","yes","y") else False
print("Low mem set to:")
print(low_mem)
L.info("Running WNN")
mu.pp.neighbors(tmp, 
                n_neighbors= int(args.n_neighbors),
                n_bandwidth_neighbors= int(args.n_bandwidth_neighbors),
                n_multineighbors= int(args.n_multineighbors),
                metric= args.metric,
                low_memory= low_mem,    
                key_added='wnn')

L.info("Computing UMAP")
mu.tl.umap(tmp,min_dist=0.4, neighbors_key='wnn')
#For taking use of mdata.obsp['connectivities'], it’s scanpy.tl.leiden() that should be used. not muon.tl.leiden
L.info("Computing Leiden clustering")
sc.tl.leiden(tmp, neighbors_key='wnn', key_added='leiden_wnn') 

L.info("Saving UMAP coordinates to csv file '%s" % args.output_csv)
umap = pd.DataFrame(tmp.obsm['X_umap'], tmp.obs.index)
umap.to_csv(args.output_csv)

#add back modalities that were not used for this run

if removed_mods is not None:
    for rmd in removed_mods:
        tmp.mod[rmd] = mdata.mod[rmd].copy()

L.info("Saving MuData to 'tmp/wnn_scaled_adata.h5mu'")
# tmp.write("tmp/wnn_scaled_adata.h5mu")
out_path = args.output_mudata
if out_path:
    if os.path.isdir(out_path) or not out_path.endswith('.h5mu'):
        out_path = os.path.join(out_path, "wnn_scaled_adata.h5mu")
else:
    out_path = "tmp/wnn_scaled_adata.h5mu"

tmp.write(out_path)

L.info("Done")

