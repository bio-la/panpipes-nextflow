#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from typing import Dict, Optional

import muon as mu
import anndata as ad

L = logging.getLogger("batch_correct_merge")
L.setLevel(logging.INFO)
_handler = logging.StreamHandler(sys.stdout)
_handler.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s - %(message)s'))
L.addHandler(_handler)



# --------------- CLI ---------------
def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument('--preprocessed_mudata', required=True, help='Base MuData (.h5mu)')
    p.add_argument('--output_mudata', required=True, help='Output MuData path (.h5mu)')

    # existing choices (keep names as your original script used)
    p.add_argument('--rna_correction_choice', default=None, help='e.g. scvi / harmony / bbknn / no_correction')
    p.add_argument('--prot_correction_choice', default=None, help='e.g. harmony / totalvi / no_correction')
    p.add_argument('--atac_correction_choice', default=None, help='e.g. harmony / bbknn / no_correction')
    p.add_argument('--multimodal_correction_choice', default=None, help='e.g. WNN / MOFA / TOTALVI')

    # NEW explicit objects (non-breaking additions)
    p.add_argument('--rna_obj', default=None, help='Path to RNA .h5ad (overrides inferred tmp/* path)')
    p.add_argument('--prot_obj', default=None, help='Path to PROT .h5ad (overrides inferred tmp/* path)')
    p.add_argument('--atac_obj', default=None, help='Path to ATAC .h5ad (overrides inferred tmp/* path)')
    p.add_argument('--multi_obj', default=None, help='Path to integrated .h5mu (WNN/MOFA/TOTALVI); overrides inferred path')
    p.add_argument('--prefer_multimodal', default='true',
                    help='true/false (default true): if a multimodal object is provided, prefer it first')

    return p


# --------------- helpers ---------------
def is_real_file(p: Optional[str]) -> bool:
    """Return True if p is a non-empty string and an existing file (and not a known sentinel)."""
    if not p:
        return False
    if os.path.basename(p) == "NO_FILE":
        return False
    return os.path.exists(p)


def merge_obs(dst: ad.AnnData, src: ad.AnnData) -> None:
    """Add new obs columns from src to dst (left-join by index). Do NOT overwrite existing columns."""
    add_cols = [c for c in src.obs.columns if c not in dst.obs.columns]
    if add_cols:
        dst.obs = dst.obs.join(src.obs[add_cols], how="left")


def merge_obsm(dst: ad.AnnData, src: ad.AnnData) -> None:
    """Copy all obsm arrays (keys present in src overwrite same keys in dst)."""
    for k, v in src.obsm.items():
        dst.obsm[k] = v


def read_unimodal(path: Optional[str]) -> Optional[ad.AnnData]:
    if is_real_file(path):
        L.info(f"Reading unimodal object: {path}")
        return ad.read_h5ad(path)
    return None


def read_multimodal(path: Optional[str]) -> Optional[mu.MuData]:
    if is_real_file(path):
        L.info(f"Reading multimodal MuData: {path}")
        md = mu.read(path)
        if not isinstance(md, mu.MuData):
            L.error(f"Provided multi_obj is not a MuData: {path}")
            sys.exit(2)
        return md
    return None


def detect_is_totalvi(mm: mu.MuData) -> bool:
    """Heuristic: treat as totalVI if any modality has obsm keys containing 'totalvi'."""
    for m in mm.mod.keys():
        if any("totalvi" in k.lower() for k in mm.mod[m].obsm.keys()):
            return True
    return False


def preserve_totalvi_obsm(mm: mu.MuData) -> Dict[str, Dict[str, object]]:
    """Collect all obsm entries whose keys contain 'totalvi' per modality from a MuData."""
    keep = {}
    for m, adata in mm.mod.items():
        keep[m] = {k: v for k, v in adata.obsm.items() if "totalvi" in k.lower()}
    return keep



parser = build_arg_parser()
args = parser.parse_args()

# base
base_object = args.preprocessed_mudata
if not is_real_file(base_object):
    L.error(f"Base MuData not found: {base_object}")
    sys.exit(2)

# Build inferred paths from choices
uni_mod_paths: Dict[str, Optional[str]] = {}

if args.rna_correction_choice and args.rna_correction_choice.lower() != "no_correction":
    uni_mod_paths['rna'] = os.path.join("tmp", f"{args.rna_correction_choice}_scaled_adata_rna.h5ad")
if args.prot_correction_choice and args.prot_correction_choice.lower() != "no_correction":
    uni_mod_paths['prot'] = os.path.join("tmp", f"{args.prot_correction_choice}_scaled_adata_prot.h5ad")
if args.atac_correction_choice and args.atac_correction_choice.lower() != "no_correction":
    uni_mod_paths['atac'] = os.path.join("tmp", f"{args.atac_correction_choice}_scaled_adata_atac.h5ad")

multi_mod = args.multimodal_correction_choice.lower() if args.multimodal_correction_choice else None
inferred_multi_path = os.path.join("tmp", f"{multi_mod}_scaled_adata.h5mu") if multi_mod else None

# Override with explicit object paths when provided
prefer_multi = str(getattr(args, 'prefer_multimodal', 'true')).lower() in ('1', 'true', 'yes', 'y')

if args.rna_obj:
    uni_mod_paths['rna'] = args.rna_obj
if args.prot_obj:
    uni_mod_paths['prot'] = args.prot_obj
if args.atac_obj:
    uni_mod_paths['atac'] = args.atac_obj

multi_path = args.multi_obj if args.multi_obj else inferred_multi_path

# Read base MuData
L.info(f"Reading base MuData: {base_object}")
base = mu.read(base_object)

present_modalities = set(base.mod.keys())
if not present_modalities:
    L.warning("Base MuData has no modalities; will still write output after merge attempts.")

# 1) Multimodal first (if present and preferred)
totalvi_keep = None
mm = None
if is_real_file(multi_path):
    mm = read_multimodal(multi_path)
    if mm is not None and prefer_multi:
        is_totalvi = (multi_mod == 'totalvi') or detect_is_totalvi(mm)
        if is_totalvi:
            L.info("Detected TOTALVI multimodal; preserving totalvi-related obsm keys before merge.")
            totalvi_keep = preserve_totalvi_obsm(mm)

        for mod in present_modalities & set(mm.mod.keys()):
            L.info(f"Merging multimodal embeddings/obs into '{mod}'")
            # First add obs columns that aren't present
            merge_obs(base.mod[mod], mm.mod[mod])
            # Then copy all obsm entries (overwrite on key collision)
            merge_obsm(base.mod[mod], mm.mod[mod])

# 2) Unimodal overlays (if provided). No prefixing (minimal change); last writer wins on obsm key collisions.
for mod, path in uni_mod_paths.items():
    if not is_real_file(path):
        continue
    if mod not in present_modalities:
        L.info(f"Base MuData missing modality '{mod}' — skipping {path}")
        continue
    L.info(f"Merging unimodal object into '{mod}': {path}")
    adata = read_unimodal(path)
    if adata is None:
        continue
    merge_obs(base.mod[mod], adata)
    merge_obsm(base.mod[mod], adata)

# If we preserved totalvi extras from mm and want to force-restore them (usually already copied above)
if totalvi_keep and mm is not None:
    L.info("Restoring preserved TOTALVI obsm keys (if they were lost).")
    for m, kv in totalvi_keep.items():
        if m in base.mod:
            for k, v in kv.items():
                base.mod[m].obsm[k] = v

base.update()
L.info(f"Writing MuData to '{args.output_mudata}'")
base.write(args.output_mudata)
L.info("Done.")
