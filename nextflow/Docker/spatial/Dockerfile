# Use Miniforge base image for ARM64 (and x86_64) compatibility
FROM condaforge/miniforge3

# Define the conda environment name
ARG conda_env=spatial_transcriptomics

# Metadata
LABEL authors="M. Gonzales A." \
      maintainer="M. Gonzales <mariana.gonzales-andre@tum.de>" \
      description="Docker image for spatial transcriptomics analysis using cell2location, scanpy, muon, spatialdata and panpipes"

# Install system dependencies (essential for compiling some Python/conda packages)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    procps \
    g++ \
    gcc \
    gfortran \
    make \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    lbzip2 \
    libgsl-dev \
    libblas-dev \
    libx11-dev \
    libreadline-dev \
    libxt-dev \
    libpcre2-dev \
    libcurl4-openssl-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy the environment YAML file into the image
COPY spatial_transcriptomics.yml /tmp/spatial_transcriptomics.yml

# Configure conda and create the environment
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels defaults && \
    conda config --set channel_priority strict && \
    conda env create -f /tmp/spatial_transcriptomics.yml && \
    conda clean -afy

# Fix for some environments requiring updated libstdc++ (especially on Apple chips)
RUN conda install -n ${conda_env} -c conda-forge libstdcxx-ng --update-deps

# Set environment variables
ENV PATH=/opt/conda/envs/${conda_env}/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/conda/envs/${conda_env}/lib
ENV CONDA_DEFAULT_ENV=${conda_env}

# Set default shell for RUN and CMD to use conda environment
SHELL ["conda", "run", "-n", "spatial_transcriptomics", "/bin/bash", "-c"]

# Run all commands inside the environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "spatial_transcriptomics"]
CMD ["/bin/bash"]
