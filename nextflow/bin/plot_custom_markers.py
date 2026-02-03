"""
Fetching data for plotting
CRG 2020-06-19
"""

import scanpy as sc
import muon as mu
from anndata import AnnData

import pandas as pd
import argparse
from itertools import product
sc.settings.autoshow = False
import os
import matplotlib
matplotlib.use('agg')
import re
from panpipes.funcs.io import  read_yaml

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
                    default="pbmc3k.h5ad",
                    help="file name, format: .h5ad")
parser.add_argument("--modalities",
                    default="rna",
                    help="modality keys containing data to be plotted, comma spearated list e.g. rna,prot default 'rna'")
parser.add_argument("--layers",
                    default="X",
                    help="layer containing data to be plotted, default=X")
parser.add_argument("--group_cols", default=["sample_id"],
                help="grouping variables", nargs='+',)
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

def relevant_group_vars(group_vars, mod):
    """
    Keep only grouping vars relevant to this modality:
        unprefixed vars apply to all modalities
        prefixed vars apply only when prefix == mod
    """
    keep = []
    for gv in group_vars:
        if ":" in gv:
            prefix, _ = gv.split(":", 1)
            if prefix == mod:
                keep.append(gv)
        else:
            keep.append(gv)
    return keep

def ensure_grouping_col(mdata, mod, gv):
    """
    Ensure mdata[mod].obs[gv] exists.
    If gv is prefixed (e.g. rna:leiden_res1), pull from mdata[mod].obs['leiden_res1'].
    Else try mdata.obs[gv], otherwise mdata[mod].obs[gv].
    """
    ad = mdata[mod]
    if gv in ad.obs.columns:
        ad.obs[gv] = ad.obs[gv].astype("category")
        return
    if ":" in gv:
        prefix, col = gv.split(":", 1)
        if prefix != mod:
            return  # not relevant to this modality
        if col not in ad.obs.columns:
            raise KeyError(
                f"Grouping var '{gv}' expects column '{col}' in mdata['{mod}'].obs, "
                f"but it was not found. Available columns: {list(ad.obs.columns)[:50]}"
            )
        ad.obs[gv] = ad.obs[col].astype("category")
        return
    # unprefixed
    if gv in mdata.obs.columns:
        ad.obs[gv] = mdata.obs.loc[ad.obs_names, gv].astype("category")
        return
    if gv in ad.obs.columns:
        ad.obs[gv] = ad.obs[gv].astype("category")
        return
    raise KeyError(
        f"Grouping var '{gv}' not found in mdata.obs or mdata['{mod}'].obs. "
        f"mdata.obs columns (head): {list(mdata.obs.columns)[:50]}; "
        f"mdata['{mod}'].obs columns (head): {list(ad.obs.columns)[:50]}"
    )

def main(adata, mod, df, grouping_var, pfx, layer_choice=None):
    for gc in df['group'].unique():
        # define file name
        # if layer_choice is None or layer_choice =="X":
        #     layer_choice = None
        #     layer_string = ""
        # else:
        #     layer_string = layer_choice
        lc = resolve_layer(adata, layer_choice, mod=mod)
        layer_string = "" if lc is None else str(lc)
        # get features
        fetches = df[df['group'] == gc]['feature']
        plot_features = [gg for gg in fetches if gg in adata.var_names]
        plot_features = list(set(plot_features))
        # do PCA if it is missing, (prerequisite for dendrogram)
        if "X_pca" not in adata.obsm.keys():
            L.warning("X_pca was not found for modality %s, but is required for the dendogram. Computing PCA with default parameters." % mod)
            sc.pp.pca(adata)
        use_dendrogram=True
        if len(adata.obs[grouping_var].unique()) > 2:
            try:
                L.info("Plotting dendogram for modality %s and group %s on X_pca" % (mod, gc))
                sc.tl.dendrogram(adata, grouping_var, use_rep="X_pca",linkage_method="average")
                use_dendrogram=True
            except ValueError: 
                use_dendrogram=False

        sc.settings.figdir  = os.path.join(args.base_figure_dir, mod ,  re.sub(":", "_", grouping_var))
        fname_prefix = "_".join([layer_string, pfx, gc ])
        fname_prefix = re.sub(":", "_", fname_prefix)
        # L.info("Plotting dotplot for modality %s, group %s, and layer %s" % (mod, gc, layer_choice))
        L.info("Plotting dotplot for modality %s, group %s, and layer %s" % (mod, gc, lc))
        sc.pl.dotplot(adata,
                        var_names=plot_features,
                        groupby=grouping_var,
                        # layer=layer_choice,
                        layer=lc,
                        dendrogram=use_dendrogram,
                        save=fname_prefix + '.png',
                        figsize=(24, 5))
        L.info("Plotting matrix plot for modality %s, group %s, and layer %s" % (mod, gc, lc))
        sc.pl.matrixplot(adata,
                        var_names=plot_features,
                        groupby=grouping_var,
                        dendrogram=use_dendrogram,
                        layer=lc,
                        save=fname_prefix + '.png',
                        figsize=(24, 5))


L.info("Reading in MuData from '%s'" % args.infile)
mdata = mu.read(args.infile)
modalities = args.modalities.split(',')

L.info("Reading in marker file from %s" % args.marker_file)
df = pd.read_csv(args.marker_file)

pfx = re.sub(".csv", "csv", os.path.basename(args.marker_file))
# write out the features that are not found in adata var (due to filtering, or incorrect name)
not_found = [gg for gg in df['feature'] if gg not in mdata.var_names]
not_found_file = re.sub(".csv", "_features_not_found.txt", os.path.basename(args.marker_file))
if len(not_found) > 0 :
    L.info("Writing features from %s that were not found in MuData to file '%s'" % (args.marker_file, not_found_file))
    with open(not_found_file, 'w') as f:
        for item in not_found:
            f.write("%s\n" % item)
            
# get layer
try:
    layers = read_yaml(args.layers)
except AttributeError:
    # this assumes that we have tried to parse a dict and nstead found a string
    # there is probably a better solution
    layers = args.layers

group_vars = args.group_cols

if type(mdata) is AnnData:
    adata = mdata
    for gv in group_vars:
        adata.obs[gv] = adata.obs[gv].astype('category')
        # main(adata, layer_choice=args.layers, group = gv, pfx = pfx, df = df)
        main(adata=adata, mod="single", df=df, grouping_var=gv, pfx=pfx, layer_choice=args.layers)
else:
    # we have multimodal object
    for mod in modalities:
        df_sub = df[df['mod'] == mod]
        # for gv in group_vars:
        #     mdata[mod].obs[gv] = mdata.obs.loc[mdata[mod].obs_names,gv].astype('category')
        # mdata.update_obs()
        # only use group vars relevant for this modalit
        group_vars_mod = relevant_group_vars(group_vars, mod)
        for gv in group_vars_mod:
            ensure_grouping_col(mdata, mod, gv)
        try:
            ll = layers[mod]
        except KeyError:
            ll = [None]
        # if len(group_vars) > 0 and ll is not None:
        #     for gv, layer in product(group_vars, ll):
        if len(group_vars_mod) > 0 and ll is not None:
            for gv, layer in product(group_vars_mod, ll):
                main(adata=mdata[mod], 
                    mod=mod,
                    layer_choice = layer,
                    grouping_var=gv, 
                    pfx = pfx,
                    df=df_sub)

L.info("Done")

