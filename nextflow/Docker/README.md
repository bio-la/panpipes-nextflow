# DOCKERFILES - Containers

This folder contains the DOCKERFILEs for the two main functionalities of this pipeline: single-cell and spatial transcriptomics.

## Overview
**Base:**  condaforge/miniforge3
**System packages:** compilers + libraries for scientific builds (gcc, g++, gfortran, zlib, lzma, bzip2, GSL, BLAS, X11, libcurl, readline, PCRE2…)
**Conda env:** single_cell created from single_cell.yml
**Conda env:** spatial_transcriptomics created from spatial_transcriptomics.yml
Fix for GLIBCXX issues → installs latest libstdcxx-ng
**ENTRYPOINT:** always runs inside the conda env (no manual activation needed)

## Build
**Single architechture :**

<pre> ```bash docker build -t repo_name/single_cell:latest . ``` </pre>

**Multi-arch (amd64 + arm64) :**
<pre> ```docker buildx create --use --name multi || true
docker buildx build --platform linux/amd64,linux/arm64 \
  -t repo_name/single_cell:latest --push . ``` </pre>

## Usage

**With Nextflow :**
<pre> ``` nextflow run main.nf -with-docker repo_name/single_cell:latest ``` </pre>