#!/usr/bin/env python
import argparse
import pandas as pd
import re
import os


from funcs.processing import  extract_parameter_from_fname

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
parser.add_argument('--input_files_str',
                    default='',
                    help='')
parser.add_argument('--output_file',
                    default='',
                    help='')
parser.add_argument('--clusters_or_markers',
                    default='',
                    help='')
parser.set_defaults(verbose=True)
args = parser.parse_args()

L.info("Running with params: %s", args)

infiles = args.input_files_str.split(" ")
# extracting sample_id
sample_ids = []
for path in infiles:
    filename = os.path.basename(path)        # e.g. "V1_Human_Lymph_Node-leiden-0.5-clusters.txt.gz"
    prefix = filename.split("-", 1)[0]      # split on first '-' only
    sample_ids.append(prefix)
sample_ids = list(set(sample_ids))

if args.clusters_or_markers == "clusters":
    L.info("Aggregating cluster columns")
    for sample in sample_ids: 
        combined_csv = pd.concat([pd.read_csv(f, sep='\t', index_col=0) for f in infiles if os.path.basename(f).startswith(sample) ], axis=1)
        # get colnames
        cnames = []
        for f in infiles:
            if os.path.basename(f).startswith(sample):
                filename = os.path.basename(f).split("-")      
                alg = filename[1]
                res = filename[2]
                cnames.append(alg + '_res_' + str(res))
                L.info(cnames)
        L.info(sample)
        L.info("Saving combined cluster columns to tsv file")
        combined_csv.to_csv("clusters/"+ sample + "-" + args.output_file, sep='\t', header=cnames, index=True)


if args.clusters_or_markers == "markers":
    L.info("Aggregating marker files")
    li = []
    all_markers_file = re.sub("_top", "_all", args.output_file)
    excel_file = re.sub("_top.txt.gz", "_all.xlsx", args.output_file)
    excel_file_top = re.sub("_top.txt.gz", "_top.xlsx", args.output_file)
    with pd.ExcelWriter(excel_file) as writer:
        with pd.ExcelWriter(excel_file_top) as writer2:
            for ff in infiles:
                # print(os.path.join(in_path, ff))
                df = pd.read_csv(ff, sep='\t')
                # add a cluster column
                clust_val = extract_parameter_from_fname(ff, "cluster", prefix=args.sample_prefix)
                sname = "cluster" + str(clust_val)
                df['cluster'] = clust_val
                li.append(df)
                df.to_excel(writer, sheet_name=sname, index=False)
                df_sub = df[df['p.adj.bonferroni'] < 0.05]
                df_sub = df_sub[df_sub['avg_logFC'] > 0]
                df_sub.to_excel(writer2, sheet_name=sname, index=False)
    frame = pd.concat(li, axis=0, ignore_index=True)
    frame.to_csv(all_markers_file, sep='\t', header=True, index=False)
    frame_sub = frame[frame['p.adj.bonferroni'] < 0.05]
    frame_sub = frame_sub[frame_sub['avg_logFC'] > 0]
    L.info("Saving combined marker files to tsv file '%s'" % args.output_file)
    frame_sub.to_csv(args.output_file, sep='\t', header=True, index=False)


L.info("Done")

