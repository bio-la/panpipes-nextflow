#!/usr/bin/env Rscript
# Clustree diagram from aggregated cluster data.
# Script originally written by Fabiola Curion, adapted by Charlotte Rich-Griffin
####


library(clustree)
library(dplyr)
library(optparse)
library(ggplot2)
library(readr)


option_list <- list(
    make_option(c("--infile"), default="none",
                help="a comma separated list of clustering files"),
	make_option(c("--plot_title"), default="none",
                help="file prefix for clustree to use as input"),
  	make_option(c("--outfile"), default="clustree.png",help = "save file")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# 1. mload and merge infiles
# 2. name columns
# 3. run clustree

message("Running with options:")

print(opt)

dir.create(dirname(opt$outfile), showWarnings = FALSE, recursive = TRUE)

# # run clustree
m = readr::read_tsv(opt$infile)
# this is a little dodge, but works ;) 
#example_column=colnames(m)[2]
#col_prefix=substr(example_column, 1, nchar(example_column)-3 )
# run clustree
#print("Running Clustree")
#gg <- clustree(m, prefix =col_prefix) + ggtitle(opt$plot_title)

prefixes <- unique(stringr::str_extract(colnames(m), "^[^_]+_(leiden|louvain)_res_"))
prefixes <- prefixes[!is.na(prefixes)]

for (prefix in prefixes) {
  print(paste("Running Clustree for:", prefix))
  gg <- clustree(m, prefix = prefix) + ggtitle(opt$plot_title)
  print(paste("Saving Clustree for:", prefix))
  ggsave(gg, filename=paste0(opt$outfile, "clustree_", prefix, ".png"), height=10,width=12, type="cairo", create.dir = TRUE)
}


#dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
#if (!(dir.exists(dirname(opt$outfile)))){
#dir.create(dirname(opt$outfile), showWarnings = FALSE, recursive = TRUE)
#}

# save
#print("Saving Clustree")
#ggsave(gg, filename=opt$outfile, height=10,width=12, type="cairo")

print("Done")
