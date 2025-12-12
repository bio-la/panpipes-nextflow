import muon as mu
import pandas as pd 
import numpy as np
import argparse
import sys
import logging
import re
import os
from muon import MuData
from anndata import AnnData


L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--input_mudata",
                    default="mdata.h5mu",
                    help="file name, format: .h5mu")
parser.add_argument("--clusters_files_csv",
                    default=None,
                    help="comma separated table of mod,compiled clusters dataframes")
parser.add_argument("--umap_files_csv",
                    default=None,
                    help = "comma separated  table of mod,txt file containing umap coordinates")
parser.add_argument("--output_mudata",
                    default="mdata.h5mu",
                    help="file name, format: .h5mu")
args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)


L.info("Reading in MuData from '%s'" % args.input_mudata)
mdata = mu.read(args.input_mudata)

 


L.info("Reading in cluster information")
cf = pd.read_csv(args.clusters_files_csv, sep='\t', index_col=0)
if args.umap_files_csv != "None": 
    L.info("Reading in UMAP coordinates")
    umap_files = args.umap_files_csv.split(" ")

if isinstance(mdata, MuData):
    pass
elif isinstance(mdata, AnnData):
    mds = cf['mod'].unique().tolist()
    if len(mds)>1:
        sys.exit("You have clustered multiple modalities but are providing only a unimodal anndata")
    else:
        L.warn("Found one modality, converting to mudata: %s " % mds[0] )    
        tmp = MuData({mds[0]:mdata})
        del mdata
        mdata = tmp
        del tmp




L.info("Adding cluster information to MuData") 
cf = cf.astype('str').astype('category')
modalities = set(col.split('_')[0] for col in cf.columns)
for modality in modalities:
    modality_cols = [col for col in cf.columns if col.startswith(modality + '_')]
    df_modality = cf[modality_cols]
    if modality != "multimodal":
        mdata[modality].obs = mdata[modality].obs.merge(df_modality, left_index=True, right_index=True)
    else:
        mdata.obs = mdata.obs.merge(df_modality, left_index=True, right_index=True)


if args.umap_files_csv != "None": 
    L.info("Adding UMAP coordinates to MuData") 
    for umap_path in umap_files: 
        modality = umap_path.split('-')[1]
        uf_df = pd.read_csv(umap_path, sep='\t', index_col=0) 
        min_dist = os.path.basename(umap_path).split("-")[2]  
        new_key = "X_umap_mindist_" + min_dist
        if modality != "multimodal": 
            if all(mdata[modality].obs_names == uf_df.index):
                mdata[modality].obsm[new_key] =  uf_df.to_numpy()
            else:
                L.warning("Cannot integrate %s into adata as obs_names mismatch" % umap_path )
        else: 
            if set(mdata.obs_names).difference(uf_df.index) == set():
                # put the observations in the same order
                uf_df = uf_df.loc[mdata.obs_names,:]
                mdata.obsm[new_key] =  uf_df.to_numpy()
            else:
                L.warning("Cannot integrate %s into mdata as obs_names mismatch" % uf.iloc[i,:] )


    
L.info("Saving updated MuData to '%s'" % args.output_mudata)
L.info(mdata)
mdata.write(args.output_mudata)
output_csv = re.sub(".h5mu", "_cell_metdata.tsv", args.output_mudata)
L.info("Saving metadata to '%s'" % output_csv)
mdata.obs.to_csv(output_csv, sep='\t')



L.info("Done")
