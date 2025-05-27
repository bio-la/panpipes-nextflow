FROM continuumio/miniconda3
ARG conda_env=preprocessing

LABEL authors="M. Gonzales A." \
  maintainer="M.Gonzales <mariana.gonzales-andre@tum.de>" \
  description="Docker image for preprocessing modules - panpipes"

# nuke cache dirs before installing pkgs
RUN rm -f /var/lib/dpkg/available && rm -rf /var/cache/apt/*

# Install system dependencies including Fortran compiler
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
    python3.11 \
    build-essential curl git procps \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Update conda
RUN conda update -n base -c defaults conda

# Add environment.yml file
COPY ingest_preprocessing.yml /tmp/ingest_preprocessing.yml

# Create conda environment and clean up
RUN conda env create -f /tmp/ingest_preprocessing.yml && \
    conda clean -afy

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "preprocessing", "/bin/bash", "-c"]

# Set PATH to use the conda environment by default
ENV PATH=/opt/conda/envs/${conda_env}/bin:$PATH
ENV CONDA_DEFAULT_ENV=${conda_env}

# Activate the environment when container starts
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "preprocessing"]
CMD ["/bin/bash"]
