from random import sample
import warnings
import re
import os
import pandas as pd
import sys
import h5py
from typing import Optional, Literal
from scanpy import read_10x_mtx, read_10x_h5, read_h5ad, read_text, read_csv, read_hdf

import muon as mu
import logging
import scirpy as ir 
from itertools import chain
from .processing import concat_adata_list, intersection
from anndata import AnnData
import numpy as np
from os import PathLike

def update_cellranger_col(path,  raw=False, method="count", sample_id=""):
    if method =="count":
        # this assumes that we get a cellranger count base path and we want file path for sc.read_10x_h5 or similar
        if raw:
            append_str = "raw_feature_bc_matrix"
        else:
            append_str = "filtered_feature_bc_matrix"
    elif method == "multi":
        if raw:
            append_str = os.path.join("multi", "count", "raw_feature_bc_matrix")
        else:
            append_str = os.path.join("per_sample_outs", sample_id, "count", "sample_filtered_feature_bc_matrix")
    else:
        sys.exit("specify cellranger filetype as cellranger_count or cellranger_multi")
    # is there a h5 file available
    if os.path.exists(os.path.join(path, append_str)) == False:
        sys.exit("cellranger path not found: %s" % os.path.join(path, append_str))
    if os.path.exists(os.path.join(path, append_str + '.h5')):
        append_str = append_str + ".h5"
        filetype = "10X_h5"
    else:
        filetype = "cellranger"
    path = os.path.join(path, append_str)
    return path, filetype


def gen_load_anndata_jobs(caf, load_raw=False, mode_dictionary = {}, load_prot_from_raw=False):
    """
    Generate a load_adatas job for each line in submission.txt
    """
    for nn in range(0, caf.shape[0]):
        if pd.isna(caf['rna_path'][nn]):
                rna_path= None
                rna_filetype=None
        elif caf['rna_filetype'][nn]=="cellranger" and mode_dictionary["rna"]:
            rna_path, rna_filetype = update_cellranger_col(caf['rna_path'][nn], raw=load_raw, method="count")
        elif caf['rna_filetype'][nn]=="cellranger_multi" and mode_dictionary["rna"]:
            rna_path, rna_filetype = update_cellranger_col(caf['rna_path'][nn], raw=load_raw, method="multi", 
                                                            sample_id=caf['sample_id'][nn])
        else:
            rna_path, rna_filetype = caf[['rna_path', "rna_filetype"]].iloc[nn]
            if load_raw:
                rna_path = re.sub("filtered", "raw", rna_path)
        # manage the prot paths
        if ('prot_path' in caf.columns and mode_dictionary["prot"]):
            # check if its the same as the rna path (data in the same file)
            if pd.isna(caf['prot_path'][nn]):
                prot_path= None
                prot_filetype=None
            elif caf['prot_path'][nn] == caf['rna_path'][nn]:
                prot_path, prot_filetype = rna_path, rna_filetype
            elif caf['prot_filetype'][nn]=="cellranger":
                # we might want to load the raw here because we want to then subset by good gex (rna) barcodes, 
                # this is why the load_prot_from_raw argument exists
                prot_path, prot_filetype = update_cellranger_col(caf['prot_path'][nn], raw=load_prot_from_raw)
            elif caf['prot_filetype'][nn]=="cellranger_multi":
                # celranger multi has the same prot and gex (rna) barcodes
                prot_path, prot_filetype = update_cellranger_col(caf['prot_path'][nn], raw=load_raw, method="multi", 
                                                                sample_id=caf['sample_id'][nn])
            else:
                prot_path, prot_filetype = caf[['prot_path', "prot_filetype"]].iloc[nn]
                if load_prot_from_raw or load_raw:
                    prot_path = re.sub("filtered", "raw", prot_path)
        else:
            prot_path= None
            prot_filetype=None
        # load tcr_path
        if 'tcr_path' in caf.columns and mode_dictionary["tcr"] and pd.notna(caf['tcr_path'][nn]):
            tcr_path = caf['tcr_path'][nn]
            tcr_filetype = caf['tcr_filetype'][nn]
        else:
            tcr_path= None
            tcr_filetype=None
        if 'bcr_path' in caf.columns and mode_dictionary["bcr"] and pd.notna(caf['bcr_path'][nn]):
            bcr_path = caf['bcr_path'][nn]
            bcr_filetype = caf['bcr_filetype'][nn]
        else:
            bcr_path= None
            bcr_filetype=None
        if ('atac_path' in caf.columns and mode_dictionary["atac"]):
            if caf.shape[0] > 1:
                sys.exit("You can only submit one atac/multiome file at a time. To aggregate, see cellranger aggr.")
            if caf['atac_filetype'][nn]=="cellranger" :
                atac_path, atac_filetype = update_cellranger_col(caf['atac_path'][nn], raw=load_raw, method="count")
            else:
                atac_path = caf['atac_path'][nn]
                atac_filetype = caf['atac_filetype'][nn]
            if 'fragments_file' in caf.columns and pd.notna(caf['fragments_file'][nn]):
                fragments_file = caf['fragments_file'][nn]
            else:
                fragments_file = None
            if 'per_barcode_metrics_file' in caf.columns and pd.notna(caf['per_barcode_metrics_file'][nn]):
                per_barcode_metrics_file = caf['per_barcode_metrics_file'][nn]
            else:
                per_barcode_metrics_file = None
            if 'peak_annotation_file' in caf.columns and pd.notna(caf['peak_annotation_file'][nn]):    
                peak_annotation_file = caf['peak_annotation_file'][nn]
            else:
                peak_annotation_file = None
        else:
            atac_path= None
            atac_filetype=None
            fragments_file = None
            peak_annotation_file = None
            per_barcode_metrics_file = None
            
        if 'barcode_mtd_path' in caf.columns:
            cell_mtd_path = caf['barcode_mtd_path'][nn]
        else:
            cell_mtd_path = None
        # create the output file 
        outfile = "./tmp/" + caf['sample_id'][nn]
        if load_raw:
            outfile = outfile + "_raw.h5mu"
        else:
            outfile = outfile + ".h5mu"
        sample_id = caf['sample_id'][nn]
        yield rna_path, outfile, \
              sample_id, \
              rna_filetype,  \
              prot_path, prot_filetype, \
              tcr_path, tcr_filetype,  \
              bcr_path, bcr_filetype, \
              atac_path, atac_filetype, \
              fragments_file, per_barcode_metrics_file, peak_annotation_file, \
              cell_mtd_path
        

def gen_load_spatial_jobs(caf, mode_dictionary = {}, load_raw=True):
    """
    Generate a load_spatial job for each line in submission.txt
    """
    for nn in range(0, caf.shape[0]):
        if "spatial_path" in caf.columns and mode_dictionary["spatial"]:
            if pd.isna(caf["spatial_path"][nn]):
                spatial_path= None
                spatial_filetype = None
            else:
                spatial_path = caf["spatial_path"][nn]
            if caf['spatial_filetype'][nn]=="xenium":
                spatial_filetype = caf['spatial_filetype'][nn]
                visium_feature_bc_matrix = None
                visium_fullres_image_file = None
                visium_tissue_positions_file = None
                visium_scalefactors_file = None
                vpt_cell_by_gene = None
                vpt_cell_metadata = None
                vpt_cell_boundaries = None
            if caf['spatial_filetype'][nn]=="vizgen":
                visium_feature_bc_matrix = None 
                visium_fullres_image_file = None
                visium_tissue_positions_file = None
                visium_scalefactors_file = None
                spatial_filetype = caf['spatial_filetype'][nn]
                vpt_cell_by_gene = None
                vpt_cell_metadata = None
                vpt_cell_boundaries = None
                if "vpt_cell_by_gene" in caf.columns:
                    if pd.notna(caf['vpt_cell_by_gene'][nn]):
                        vpt_cell_by_gene = caf['vpt_cell_by_gene'][nn]
                if "vpt_cell_metadata" in caf.columns:
                    if pd.notna(caf['vpt_cell_metadata'][nn]):
                        vpt_cell_metadata = caf['vpt_cell_metadata'][nn]
                if "vpt_cell_boundaries" in caf.columns:
                    if pd.notna(caf['vpt_cell_boundaries'][nn]):
                        vpt_cell_boundaries = caf['vpt_cell_boundaries'][nn]
            elif caf['spatial_filetype'][nn]=="visium":
                visium_feature_bc_matrix = None 
                visium_fullres_image_file = None
                visium_tissue_positions_file = None
                visium_scalefactors_file = None
                vpt_cell_by_gene = None
                vpt_cell_metadata = None
                vpt_cell_boundaries = None
                spatial_filetype = caf['spatial_filetype'][nn]
                #counts file
                if "visium_feature_bc_matrix" in caf.columns:
                    if pd.notna(caf["visium_feature_bc_matrix"][nn]):
                        visium_feature_bc_matrix= caf["visium_feature_bc_matrix"][nn]
                # fullres image
                if "visium_fullres_image_file" in caf.columns:
                    if pd.notna(caf["visium_fullres_image_file"][nn]):
                        visium_fullres_image_file= caf["visium_fullres_image_file"][nn]
                # tissue position 
                if "visium_tissue_positions_file" in caf.columns:
                    if pd.notna(caf["visium_tissue_positions_file"][nn]):
                        visium_tissue_positions_file= caf["visium_tissue_positions_file"][nn]
                # scalefactor
                if "visium_scalefactors_file" in caf.columns:
                    if pd.notna(caf["visium_scalefactors_file"][nn]):
                        visium_scalefactors_file= caf["visium_scalefactors_file"][nn] 
        else:
            spatial_path= None
            spatial_filetype = None
            visium_feature_bc_matrix = None
            visium_fullres_image_file = None
            visium_tissue_positions_file = None
            visium_scalefactors_file = None
            vpt_cell_by_gene = None
            vpt_cell_metadata = None
            vpt_cell_boundaries = None
            
        if 'barcode_mtd_path' in caf.columns:
            cell_mtd_path = caf['barcode_mtd_path'][nn] #not yielding this right now!
        else:
            cell_mtd_path = None
        # create the output file 
        outfile = "./tmp/" + caf['sample_id'][nn]
        if load_raw:
            outfile = outfile + "_raw.zarr" 
        else:
            outfile = outfile + ".zarr"
        sample_id = caf['sample_id'][nn]

        yield spatial_path, outfile, sample_id, spatial_filetype, \
              visium_feature_bc_matrix, visium_fullres_image_file, visium_tissue_positions_file, visium_scalefactors_file, \
              vpt_cell_by_gene, vpt_cell_metadata, vpt_cell_boundaries


def read_anndata(
    fname: Optional[str] = None,
    use_muon: Optional[bool] = False, 
    modality: Literal['all', 'rna', 'prot','atac','rep'] = 'all'):
    logging.info("reading %s" % fname)
    # check fname file exists
    try:
        os.path.isfile(fname)
    except FileNotFoundError:
        sys.exit("anndata file not found")
    if use_muon is False:
        return read_h5ad(fname)
    else:
        logging.info("reading %s modality" % modality)
        if modality=="all":
            # note this just loads the rna part of the muon object
            return mu.read(fname)
        if modality=="rna":
            # note this just loads the rna part of the muon object
            return mu.read(fname + "/rna")
        elif modality=="prot":
            return mu.read(fname + "/prot")
        elif modality=="atac":
            return mu.read(fname + "/atac")
        elif modality=="rep":
            return mu.read(fname + "/rep")
        
        else:
            sys.exit("modality not found, must be one of 'all', 'rna', 'prot','atac', 'rep' ")


def write_anndata(adata,
                    fname, 
                    use_muon=False, 
                    modality: Literal['all', 'rna', 'prot'] = 'all'):
    if use_muon is False:
        adata.write(fname)
    else:
        if modality=="rna":
            # note this just saves the rna part of the muon object
            return mu.write(fname + "/rna", adata)
        if modality=="prot":
            # note this just save the prot part of the muon object
            return mu.write(fname + "/prot", adata)
        elif modality=="all":
            return mu.write(fname, adata)
        else:
            sys.exit("modality not found, must be one of 'all', 'rna', 'prot'")

# --- loading raw dat

def write_obs(mdata, output_prefix="./", output_suffix="_filtered_cell_metadata.tsv"):
        metafile = mdata.obs
        metafile["cellbarcode"] = mdata.obs.index
        savename= output_prefix + output_suffix
        metafile.to_csv(savename, sep='\t', index=True)

def check_submission_file(caf):
    if "filetype" in caf.columns or "path" in caf.columns:
        raise ValueError("you appear to be using the old notation for the sample submission file, \
        please update to use rna_path instead of path and rna_filetype instead of filetype")
    # check for required cols
    req_cols = ['sample_id', 'rna_path', 'rna_filetype']
    for rc in req_cols:
        if rc not in caf.columns:
            raise ValueError("required column %s missing from submission file" % rc)



def check_filetype(path, filetype):
    logging.debug(path)
    logging.debug(filetype)
    if filetype == "cellranger" :
        if os.path.exists(path + "/matrix.mtx.gz"):
            logging.info("cellranger matrix found")
        elif os.path.exists(path + "/matrix.mtx"):
            logging.info("cellranger matrix found")
        else:
            sys.exit("cellranger matrix file not found, have you specified the path correctly")
    else:
        ftype_checks = {
            "h5ad":r".h5ad$",
            "csv_matrix": ".csv",
            "txt_matrix" : ".txt",
            "10X_h5": ".h5",
            "hdf": ".h5",
            "cellranger_vdj": r".json|.csv",
            "vizgen": r".txt|.csv|.tsv" #suboptimal for now but roll with it
        }
        if filetype not in ftype_checks.keys():
            sys.exit("unknown filetype %s, please specify one of: %s" % filetype, ftype_checks.keys().join(", "))
        elif re.search(ftype_checks[filetype], path) is None:
            sys.exit("filetype does not match file suffix, check it is specified correctly \
                    (with %s as suffix (gzipped files ok for csv and tsv))" % ftype_checks[filetype])
        elif os.path.exists(path)==False:
            sys.exit("file %s does not exist" % path)
        else:
            logging.info("%s file exists" % path)


def read_scirpy(
    fname, 
    filetype: Literal['cellranger_vdj', 'tracer', 'bracer', 'airr' ] = 'cellranger_vdj'):
    """
    wrapper around all the scirpy io functions
    """
    # check fname file exists
    logging.info("loading repertoire using scirpy, filetype %s" % filetype)
    try:
        os.path.isfile(fname)
    except FileNotFoundError:
        sys.exit("vdj file not found: %s" % fname)
    if filetype ==  "cellranger_vdj":
        vdj = ir.io.read_10x_vdj(fname)
    elif filetype == "tracer":
        vdj = ir.io.read_tracer(fname)
    elif filetype ==  "bracer":
        vdj = ir.io.read_bracer(fname)
    elif filetype ==  "airr":
        vdj = ir.io.read_airr(fname)
    else:
        sys.exit("fileype not valid, must be one of 'cellranger_vdj', 'tracer', 'bracer', 'airr'")
    return vdj


def scp_read_10x_h5(filename: PathLike, library_keep=None, extended= False, *args, **kwargs) -> AnnData:
    """
    expanded sc.read_10x_h5 to filter for the library 
    adapted from https://github.com/scverse/muon/blob/master/muon/_prot/io.py
    """
    adata = read_10x_h5(filename, gex_only=False, *args, **kwargs)
    logging.debug(adata)
    logging.debug("filtering cellranger outputs to %s" % library_keep)
    if extended is not None and extended is True:
        import shutil
        logging.debug("copying %s" % filename )
        h5file = h5py.File(filename, "r")
        logging.debug("reading intervals")
        if "interval" in h5file["matrix"]["features"]:
            intervals = np.array(h5file["matrix"]["features"]["interval"]).astype(str)
            logging.debug("intervals %s" %intervals)
            adata.var["interval"] = intervals
            logging.debug(f"Added `interval` annotation for features from {filename}")
            logging.debug("updated adata %s" % adata)
            logging.debug("adata.var is %s" % adata.var)
            h5file.close()
        else:
            # Make sure the file is closed
            h5file.close()
    else:
        logging.debug("skipping extended arg")
    if library_keep is not None:
        adata = adata[:, list(map(lambda x: x == library_keep, adata.var["feature_types"]))].copy()
        logging.debug("filtered anndata size %s" % adata)
    
    return adata



def scp_read_10x_mtx(filename: PathLike, library_keep=None, *args, **kwargs) -> AnnData:
    """
    expanded sc.read_10x_mtx to filter for the library 
    adapted from https://github.com/scverse/muon/blob/master/muon/_prot/io.py
    """
    adata = read_10x_mtx(filename, gex_only=False, *args, **kwargs) #need to leave gex as this is scanpy's code
    logging.debug("filtering cellranger outputs to %s" % library_keep)
    if library_keep is not None:
        adata = adata[
            :, list(map(lambda x: x == library_keep, adata.var["feature_types"]))
        ].copy()
    return adata

#----------------
# temp loading functions for spatial data
# load_Spatial_in
# load_spatial_from_multiple_files
# ->>>how to write?
#----------------





def load_adata_in(path, filetype, gex_only=True, var_names="gene_symbols", library=None, extended=False):
    """
    load in any format supported by the load_functions_dict
    expected to be called by load_mdata_from_multiple_files
    """
    # all the supported loading functions here
    load_functions_dict = {
    "cellranger": lambda x: scp_read_10x_mtx(x, 
                                           library_keep=library,
                                           var_names=var_names,
                                           cache=False),
     "h5ad":read_h5ad,
     "csv_matrix":read_csv,
     "txt_matrix": read_text,
     "10X_h5": lambda x: scp_read_10x_h5(x, library_keep=library, extended=extended),
     "hd5": read_hdf,
      # vdj functions
     "cellranger_vdj": lambda x: read_scirpy(x, filetype="cellranger_vdj"),
     "tracer":lambda x: read_scirpy(x, filetype="tracer"),
     "bracer":lambda x: read_scirpy(x, filetype="bracer"),
     "airr":lambda x: read_scirpy(x, filetype="airr"),             
     }
    

    # the var_names argument decides which column of the cellranger features files ends up as the index.
    if var_names=="gene_symbols":
        col_fill="gene_ids"
    elif var_names=="gene_ids":
        col_fill="gene_symbols"
    else:
        col_fill="index"

    # try:
    logging.debug("loading %s " % filetype)
    adata = load_functions_dict[filetype](path)
    logging.debug("this is the anndata now: %s" % adata)
    # in some cases you need to update the var index col.
    if var_names in adata.var.columns:
        logging.info("resetting var.index")
        # we need to reset the index to be the expected one.
        adata.var = adata.var.reset_index().rename(columns={"index":col_fill}).set_index(var_names)
    return adata
    # except:
    #     raise ValueError("unknown filetype %s, please specify one of: %s" % filetype, ", ".join(list(load_functions_dict.keys()))) 
    


def update_intersecting_feature_names(adata1, adata2, prefix):
    # check for intersecting names
    if (len(adata1.var_names) > 0 and len(adata2.var_names) > 0 ):
        intersect_features = intersection(adata1.var_names.tolist(), adata2.var_names.tolist())
        if len(intersect_features) > 0:
            # this will be a problem downstream, so we prefix the 
            warnings.warn("intersecting gene ids and other modalities are problematic downstream, \
                prefixing  with %s" % prefix)
            adata2.var[ prefix + 'id'] = adata2.var_names
            adata2.var_names = [prefix + x for x in adata2.var_names.tolist()]
    return adata2


def merge_tcr_bcr_into_one_anndata(tcr, bcr):
    logging.info("merging tcr and bcr into one rep modality")
    intersect_obs = list(set(bcr.obs_names).union(set(tcr.obs_names)))
    adata = AnnData(X=np.empty(shape=(len(intersect_obs), 0)),
                       obs=pd.DataFrame(index=intersect_obs))
    # merge in the tcr and bcr with the rna 
    ir.pp.merge_airr_chains(adata , tcr)
    ir.pp.merge_airr_chains(adata , bcr)
    ir.tl.chain_qc(adata)
    return adata
    


def _make_one_rep_modality(data_dict):
    """
    takes the data dict from load_mdata_from_multiple_files, and collapses down the bcr and/or tcr into one modality called rep
    """
    if ("tcr" in data_dict.keys() and "bcr" in data_dict.keys()):
        # get the intersecting obs
        data_dict['rep'] = merge_tcr_bcr_into_one_anndata(data_dict['tcr'], data_dict['bcr'])
        # remove separate assays
        data_dict.pop('tcr', None)
        data_dict.pop('bcr', None)
    elif("tcr" in data_dict.keys()):
        logging.info("storing tcr in rep modality")
        data_dict['rep'] = data_dict['tcr']
        data_dict.pop('tcr', None)
    elif("bcr" in data_dict.keys()):
        logging.info("storing bcr in rep modality")
        data_dict['rep'] = data_dict['bcr']
        data_dict.pop('bcr', None)
    


def load_mdata_from_multiple_files(all_files_dict):
    """
    create a mudata object from multiples files
    Parameters
    ----------
    all_files_dict: dict
        dictionary containing one key per assay out of [RNA, PROT, TCR, BCR, ATAC]
        and the values for each key is a list of file path and fie type
        e.g.  {"RNA": [filepath, "filetype"],
               "PROT: [file path, "filetype"]} 
        Filetypes supported for RNA/prot: ["cellranger", "h5ad", "csv_matrix", "txt_matrix", "10X_h5"],
        Filetypes supported for atac (multiome preferred is 10X_h5) ["10X_h5","cellranger","h5ad"]
        Filetypes supported for rep: ["cellranger_vdj", "airr", "tracer", "bracer"  ] 
        See scirpy documentation for more information of repertoire input formats 
        https://scverse.org/scirpy/latest/api.html#module-scirpy.io
    """
    # convert names to match mudata conventions
    # mudata_conventional_names
    # {"GEX":"rna", "ADT":"prot", "TCR":"tcr", "BCR":"bcr", "ATAC": "atac", "spatial":"spatial"}
    # all_files_dict = {mudata_conventional_names[nm]: x  for (nm, x) in all_files_dict.items()}
    # note: scanpy's default function use gex_only as param so we need to leave that in
    logging.debug(all_files_dict.keys())
    # load in separate anndata for each expected modality
    data_dict = {}
    for nm, x in all_files_dict.items():
        logging.info("loading %s" % nm)
        extra_args={}
        if nm == "rna":
            extra_args['library'] = "Gene Expression"
            if 'atac' in all_files_dict.keys():
                extra_args["extended"]= True
            else:
                extra_args["extended"]= False
        if nm == "prot":
            extra_args["gex_only"] = False
            extra_args["var_names"] = "gene_ids"
            extra_args['library'] = "Antibody Capture"
            extra_args["extended"]= False
        if nm == "atac":
            extra_args["gex_only"] = False
            extra_args['library'] = "Peaks"
            extra_args['extended'] = False
        if nm == "spatial":
            extra_args["gex_only"] = True # check this for techs other than merfish and visium H&E
            #extra_args["counts_file"] =
            extra_args['extended'] = False

        logging.debug("extra args")
        logging.debug(extra_args)
        data_dict[nm] = load_adata_in(x[0], x[1], **extra_args) #** 
        #x[0] is the path, x[1] is the filetype
    logging.debug(data_dict["rna"])
    logging.debug(data_dict.keys())
    # we want unique var names for each assay
    if "rna" in data_dict.keys():
        logging.debug("rna adata %s " % data_dict["rna"])
        data_dict['rna'].var_names_make_unique()
        for mod in data_dict.keys():
            if mod =="rna":
                # skip since we want to update all the ones that are not rna
                pass
            else:
                logging.debug("second adata %s " % data_dict[mod])
                logging.debug("updating intersection for %s" % mod)
                data_dict[mod] = update_intersecting_feature_names(adata1=data_dict['rna'], adata2=data_dict[mod], 
                prefix=mod + "_")
    #make repertoire data into one assay
    _make_one_rep_modality(data_dict)
    logging.debug(data_dict)

    # finally make the mudata object
    mdata_ = mu.MuData(data_dict)
    mdata_.update()
    # remove empty_rna
    return mdata_


from scipy.io import mmwrite
import gzip
def write_10x_counts(adata, path, layer=None):
    if layer is None:
        arr = adata.X
    else:
        arr = adata.layers[layer]
    if os.path.exists(path) is False:
        os.makedirs(path)
    # convert to integer if appropriate
    if arr.sum().is_integer():
        arr = arr.astype(int)
    # write matrix.mtx
    mmwrite(os.path.join(path, "matrix.mtx"), arr.T)
    with open(os.path.join(path, "matrix.mtx"), 'rb') as f_in, gzip.open(os.path.join(path, "matrix.mtx.gz"), 'wb') as f_out:
        f_out.writelines(f_in)
    os.remove(os.path.join(path, "matrix.mtx"))
    # write barcodes file
    barcodes = adata.obs_names.to_frame()
    barcodes.to_csv(os.path.join(path, "barcodes.tsv.gz"), sep='\t', index=None, header=None)
    # write features file (trying to make sure the correct columns are stored)
    if "gene_ids" in adata.var.columns:
        # we assume the index is gene_ids
        features = adata.var.reset_index().rename(columns={"index": "gene_symbols"})
    elif "gene_symbols" in adata.var.columns:
        features = adata.var.reset_index().rename(columns={"index": "gene_ids"})
        if "prot_id" in features.columns:
            features = features.drop(columns=['gene_symbols']).rename(columns={"prot_id" : "gene_symbols"})
    features = features[['gene_ids','gene_symbols', 'feature_types']] 
    features.to_csv(os.path.join(path, "features.tsv.gz"), sep='\t', index=None, header=None)

# https://stackoverflow.com/questions/33529312/remove-empty-dicts-in-nested-dictionary-with-recursive-function
def dictionary_stripper(data):
    new_data = {}
    # Only iterate if the given dict is not None
    if data:
        for k, v in data.items():
            if isinstance(v, dict):
                v = dictionary_stripper(v)
            # ideally it should be not in, second you can also add a empty list if required
            if v not in ("", None, {}, []):
                new_data[k] = v
        # Only if you want the root dict to be None if empty
        if new_data == {}:
            return None
        return new_data
    return None


def replace_string_nones(d):
    for k, v in d.copy().items():
        if isinstance(v, dict):
            replace_string_nones(v)
        elif v == "None":
            d[k] = None

import yaml
def read_yaml(fname):
    if type(fname) is dict:
        params=fname
    elif os.path.exists(fname):
        with open(fname, "r") as stream:
            try:
                params = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
    else:
        # if the yaml is being passed as a string
        try:
            params = yaml.safe_load(fname)
            
        except yaml.YAMLError as exc:
            print(exc)
    replace_string_nones(params)
    return params


# def model_load_choice( ):
#     import scvi
#     functions_dict={
#     "totalvi": lambda x: scvi.model.TOTALVI.load_query_data,
#     "scvi": lambda x: scvi.model.SCVI.load_query_data,
#     "scanvi": lambda x: scvi.model.SCANVI.load_query_data,

#     }
    
