"""
Plotting umaps of a small list of custom genes
CRG 2020-06-19
"""
import scanpy as sc
import muon as mu
from anndata import AnnData
import pandas as pd
import argparse
from panpipes.funcs.io import read_yaml
from itertools import product

sc.settings.autoshow = False
import os
import matplotlib
matplotlib.use('agg')

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--infile",
                    default="pbmc3k.h5mu",
                    help="file name, format: .h5mu/.h5ad")
parser.add_argument("--modalities",
                    default="rna",
                    help="modality keys containing data to be plotted, comma spearated list e.g. rna,prot default 'rna'")
parser.add_argument("--layers",
                    default="X",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--basis_dict",
                    default="X_umap",
                    help="basis in obsm to plot")
parser.add_argument("--marker_file",
                    default="Marker_Lists-myeloid.csv",
                    help="comma separateed list of marker_list csvs, one column of gene_ids")
parser.add_argument("--base_figure_dir", default="figures/",
                    help="figures path")

args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)
sc.settings.figdir = args.base_figure_dir

# ---- script
def resolve_layer(adata, layer_choice, mod=None):
    """
    Return a valid layer key for this AnnData, or None to use .X.
        Treat None / "" / "X" as "use .X"
        If the object has no layers at all, always use .X
        If a layer is requested but missing, warn and fall back to .X
    """
    if layer_choice in (None, "", "X"):
        return None
    # If there are no layers stored, must use X
    if len(getattr(adata, "layers", {}).keys()) == 0:
        return None
    if layer_choice not in adata.layers.keys():
        L.warning(
            "Layer '%s' not found in modality %s. Available layers: %s. Falling back to .X",
            layer_choice,
            mod if mod is not None else "<unknown>",
            list(adata.layers.keys())
        )
        return None
    return layer_choice

def main(adata, mod, layer_choice, df, basis):
    for gc in df['group'].unique():
        # define file name
        # if layer_choice is None or layer_choice =="X":
        #     layer_choice = None
        #     layer_string = ""
        # else:
        #     layer_string = layer_choice
        # fname_prefix = "_".join(["_" + mod, layer_string, gc])
        # resolve layer safely for this modality
        # use the helper function to get the correct layer
        lc = resolve_layer(adata, layer_choice, mod=mod)
        layer_string = "" if lc is None else str(lc)
        fname_prefix = "_".join(["_" + mod, layer_string, gc]).strip("_")
        # get features
        fetches = df[df['group'] == gc]['feature']
        plot_features = [gg for gg in fetches if gg in adata.var_names]
        plot_features = list(set(plot_features))
        L.info("Plotting embedding %s of modality %s with layer %s" % (basis, mod, layer_choice))
        sc.settings.figdir  = os.path.join(args.base_figure_dir, mod)
        #mu.pl.embedding(adata, basis=basis, layer=layer_choice, color=plot_features, save = fname_prefix + ".png")
        mu.pl.embedding(adata, basis=basis, layer=lc, color=plot_features, save=fname_prefix + ".png")

L.info("Reading in MuData from '%s'" % args.infile)
mdata = mu.read(args.infile)
modalities = args.modalities.split(',')

L.info("Reading in marker file from %s" % args.marker_file)
df = pd.read_csv(args.marker_file )

# get layer
try:
    layers = read_yaml(args.layers)
except AttributeError:
    # this assumes that we have tried to parse a dict and nstead found a string
    # there is probably a better solution
    layers = args.layers

# get bases
try:
    basis_dict = read_yaml(args.basis_dict)
except AttributeError:
    # this assumes that we have tried to parse a dict and nstead found a string
    # there is probably a better solution
    basis_dict = args.basis_dict


if type(mdata) is AnnData:
    adata = mdata
    #basis= basis_dict
    #main(adata, layer_choice=args.layers,  df = df, basis=basis)
    # basis_dict may be a string here (single basis)
    basis = basis_dict if isinstance(basis_dict, str) else "X_umap"
    main(adata=adata, mod="single", layer_choice=args.layers, df=df, basis=basis)
else:
    # we have multimodal object
    for mod in modalities:
        print(mod)
        df_sub = df[df['mod'] == mod]
        mdata.update_obs()
        try:
            ll = layers[mod]
        except KeyError:
            ll = [None]
        # if mod in basis_dict.keys():
        #     bb= basis_dict[mod]
        # else:
        #     bb = []
        # basis_dict can be either a dict per modality, or a single string basis
        if isinstance(basis_dict, dict):
            bb = basis_dict.get(mod, [])
        else:
            bb = [basis_dict]
        if len(bb) > 0 :
            for basis, layer in product(bb, ll):
                print(basis,layer)
                if basis.startswith("X_"):
                    key = basis
                else:
                    key = f"X_{basis}"
                available = mdata[mod].obsm_keys()
                if key not in available and basis not in available:
                    L.warning("Basis '%s' not found in mdata['%s'].obsm; skipping. Available: %s",basis, mod, list(available))
                    continue
                main(adata=mdata[mod], 
                    mod=mod,
                    layer_choice = layer,
                    df=df_sub,
                    basis=basis)

L.info("Done")

