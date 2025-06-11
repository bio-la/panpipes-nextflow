# Use Miniforge base image for ARM64 (and x86_64) compatibility
FROM condaforge/miniforge3

ARG conda_env=preprocessing

LABEL authors="M. Gonzales A." \
      maintainer="M.Gonzales <mariana.gonzales-andre@tum.de>" \
      description="Docker image for preprocessing modules - panpipes"


# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential curl git procps \
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy conda environment YAML
COPY ingest_preprocessing.yml /tmp/ingest_preprocessing.yml
# Add conda channels and set strict priority
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels defaults && \
    conda config --set channel_priority strict

# Create conda environment and clean up
RUN conda env create -f /tmp/ingest_preprocessing.yml && \
    conda clean -afy

# Ensure latest libstdcxx-ng in conda env for GLIBCXX compatibility
RUN conda install -n ${conda_env} -c conda-forge libstdcxx-ng --update-deps

# Set environment variables for conda environment
ENV PATH=/opt/conda/envs/${conda_env}/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/conda/envs/${conda_env}/lib:$LD_LIBRARY_PATH
ENV CONDA_DEFAULT_ENV=${conda_env}

# Use conda environment by default in RUN commands
SHELL ["conda", "run", "-n", "preprocessing", "/bin/bash", "-c"]

# Activate environment on container start
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "preprocessing"]
CMD ["/bin/bash"]
