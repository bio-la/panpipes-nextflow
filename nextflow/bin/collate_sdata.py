#!/usr/bin/env python
import pandas as pd 
import argparse
import sys
import logging
import spatialdata as sd
import os


L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

# parse argumetns
parser = argparse.ArgumentParser()
parser.add_argument("--input_sdata",
                    help="List of sdatas")
parser.add_argument("--clusters_files_csv",
                    default=None,
                    help="List of aggregated cluster files")
parser.add_argument("--umap_files_csv",
                    default=None,
                    help = "List of UMAP coords files")
parser.add_argument("--output_dir",
                    default="./clustered.data/",
                    help="Output directory")
args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)


sample_paths = args.input_sdata.split(" ")
samples_dic = {}
for path in sample_paths: 
    filename = os.path.basename(path)        
    prefix = filename.split("-", 1)[0]     
    samples_dic[prefix] = sd.read_zarr(path)

cluster_files = args.clusters_files_csv.split(" ")
cluster_dic = {}
for path in cluster_files: 
    filename = os.path.basename(path)        
    prefix = filename.split("-", 1)[0]     
    cluster_dic[prefix] = path

umap_files = args.umap_files_csv.split(" ")


# add in the clusters
for sample_id, sdata in samples_dic.items():
    L.info("Adding Cluster information to the SpatialData of sample %s", sample_id)
    cluster_file = cluster_dic[sample_id]
    L.info(cluster_file)
    cf_df = pd.read_csv(cluster_file, sep='\t', index_col=0) 
    cf_df = cf_df.astype('str').astype('category')
    sdata["table"].obs = sdata["table"].obs.merge(cf_df, left_index=True, right_index=True)
    samples_dic[sample_id] = sdata



L.info("Adding UMAP coordinates to SpatialData")

for sample_id, sdata in samples_dic.items():
    L.info("Adding UMAP coordinates to the SpatialData of sample %s", sample_id)
    for umap_path in umap_files: 
        if umap_path.startswith(sample_id):
            L.info(umap_path)
            uf_df = pd.read_csv(umap_path, sep='\t', index_col=0) 
            min_dist = os.path.basename(umap_path).split("-")[1]   
            if all(sdata["table"].obs_names == uf_df.index):
                new_key = "X_umap_mindist_" + min_dist
                sdata["table"].obsm[new_key] =  uf_df.to_numpy()
            else:
                L.warning("Cannot integrate %s into adata as obs_names mismatch" % umap_path )
    samples_dic[sample_id] = sdata
    



L.info("Saving updated SpatialData to '%s'" % args.output_dir)
for sample_id, sdata in samples_dic.items():
    sdata.write(args.output_dir + sample_id + ".zarr")
    L.info("Saving metadata to " + sample_id + "_cell_metdata.tsv")
    sdata["table"].obs.to_csv(sample_id + "_cell_metdata.tsv", sep='\t')


L.info("Done")
