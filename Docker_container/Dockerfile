FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

LABEL maintainer="aurelio.maiolino@edu.unito.it"

# Install 
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gpg \
    curl \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    make \
    cmake \
    gfortran \
    libxt-dev \
    liblapack-dev \
    libblas-dev \
    sudo \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    pandoc \
    git

# Cleaning step
RUN rm -rf /var/lib/apt/lists/*



# Add the CRAN GPG key and repository for R
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    | tee /etc/apt/sources.list.d/cran-r.list

# Update again and install R
RUN apt update && apt install -y --no-install-recommends r-base


# Install JupyterLab
RUN apt update && apt install -y python3 python3-pip python3-venv

# create a virtual environment in which JupyterLab can be installed
RUN python3 -m venv /opt/venv

# Activate virtual environment and install JupyterLab
RUN /opt/venv/bin/pip install --upgrade pip && /opt/venv/bin/pip install jupyterlab

# Set the virtual environment as the default Python path
ENV PATH="/opt/venv/bin:$PATH"

# Make R visible to jupyter
RUN R -e "install.packages('IRkernel')" \
    R -e "IRkernel::installspec(user = FALSE)"


# Install R packages
RUN R -e "install.packages(c('BiocManager', 'dplyr', 'ggplot2', 'data.table', 'future', 'cowplot', 'remotes', 'R.utils', 'dplyr', 'rtracklayer', 'tinytex'))" 
RUN R -e "BiocManager::install('tidyverse')" 

# Install Seurat
RUN R -e "BiocManager::install('Seurat')"

# Install Signac
RUN R -e "remotes::install_github('stuart-lab/signac', ref = 'develop')"

# Install custom R packages
RUN R -e "remotes::install_github('Maiolino-Au/PABMaiolinoPackage')"

# Download the R script from GitHub
RUN mkdir -p /0_Scripts && \
  cd /0_Scripts && \
  curl -O https://raw.githubusercontent.com/Maiolino-Au/ProgrAppBioinfo_Exam/main/Scripts/Script_PAB_exam_Maiolino.R && \
  curl -O https://raw.githubusercontent.com/Maiolino-Au/ProgrAppBioinfo_Exam/main/Scripts/Maiolino_Au.Rmd

# 
ENV SHELL=/bin/bash
CMD jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root --ServerApp.allow_origin='*' --ServerApp.token='' #last one disables the token
