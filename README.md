# Programming Approaches for Bioinformatics Exam
Maiolino Aurelio - June 2025 \
Exam of the course Programming Approaches for Bioinformatics 

## Index
1. [Project](#project)
2. [Docker](#docker)
3. [Logic](#logic)
4. [other](#other)




## Project 
Deadline: 26/06/2025

### Part 1: Docker Setup
* Create a Docker container using a Dockerfile that starts from a base Ubuntu image.

### Part 2: Package Installation
* Attempt to install the following R packages:
   * Seurat (follow [installation guide](https://satijalab.org/seurat/articles/install_v5.html))
   * Signac ([GitHub repository](https://github.com/stuart-lab/signac))
* If any installation fails, in the report, provide a clear explanation of the issues encountered and why the installation was unsuccessful.
* Possibly try to find a work around to the installation fails, if you fail at all in installing the package(s) describe the reason why it was impossible to find a work around.

### Part 3: Data Processing and Analysis
Use the dataset: "PBMCs 3k cells from a healthy donor"
* Use the [material](https://drive.google.com/file/d/1Yc5kv8OXd5oDTkX-TVRvwLkwkjew43DZ/view?usp=drive_link) provided as part of the exam.

#### Step 1: Matrix Conversion
* Use the Matrix R package to convert the sparse matrix into a full matrix.
* Save the result as a data.table object.

#### Step 2: Split Gene Expression and ATAC-seq Data
From the data.table object, separate:
* Gene expression data (rows labelled with Ensembl gene IDs, e.g., ENSG00000243485)
* ATAC-seq peak data (rows labelled with genomic coordinates, e.g., chrN:NNNN-NNNN)

#### Step 3: Summarize Data
For each dataset (expression and peaks), compute the column-wise sum to produce:
* A single vector of total expression per gene
* A single vector of total chromatin accessibility per peak region

#### Step 4: Create Genomic Ranges
* Convert both the summarized gene expression and peak data into GenomicRanges objects.
* Add the summarized data as metadata to their respective GenomicRanges.

#### Step 5: Gene Annotation for ATACseq data
* Using the annotation file Homo_sapiens.GRCh38.114.gtf.gz:
   * Create a GenomicRanges object only for protein-coding genes and only for gene features.
   * Remap the ATAC-seq GenomicRanges to this object and attach the summarized peak data from step 4.

#### Step 6: Finalize Expression Data
* Subset the expression GenomicRanges, step 4, to include only protein-coding genes.
* Add gene symbol identifiers to the object.

#### Step 7: Data Normalization and Integration
* Normalize both expression and ATAC-seq data using CPM:
   * Divide each column by the column sum, multiply by 106, add a pseudo-count of 1, and apply log2.
* Merge expression and ATAC data based on common genes.
* Provide a summary table of the number of ATAC peaks that could not be merged and a plot of peak intensity distribution chromosome by chromosome. Provide a summary table of the genes which do not show association with ATAC peaks and plot their expression distribution chromosome by chromosome

#### Step 8: Visualization
* Generate a scatter plot using ggplot2:
   * X-axis: log-transformed expression CPM
   * Y-axis: log-transformed ATAC CPM
* If the plot is too much busy of data divide the plot in the 24 chromosomes

### Part 4: Reporting
* Create an HTML report using RMarkdown that includes:
   * All code used
   * Output results (tables, plots, etc.)
   * A description of any problems or errors encountered
   * Possible workarounds or solutions you tried
   * Describe how to repeat the analysis using the docker provided.
* Provide the docker with the software installed which includes the scripts used to perform the analysis.




## Docker
I used a docker container with stuff

### Dockerfile
aaa

Build the container from satijalab/seurat v5.0.0. This is an image which already contains R and Seurat, it is quicker compared to downloading them singularly and enshure the use of spcecific versions of the various libraries. \
R v4.2.0 \
Seurat v5.0.0
```Dockerfile
FROM ubuntu
```

Setup
```Dockerfile
ENV DEBIAN_FRONTEND=noninteractive

LABEL maintainer="aurelio.maiolino@edu.unito.it"
```


Install curl: it is needed to download the automatically scripts. last step is cleaning
```Dockerfile
# Install stuf for ubuntu
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
    && rm -rf /var/lib/apt/lists/*
```

Intstall R
```Dockerfile
# Add the CRAN GPG key and repository for R
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    | tee /etc/apt/sources.list.d/cran-r.list

# Update again and install R
RUN apt update && apt install -y --no-install-recommends r-base
```

Install Python, Jupyter, and make R visible by jupiter.
```Dockerfile
# Install JupyterLab
RUN apt update && apt install -y python3 python3-pip python3-venv \
    python3 -m venv /opt/venv \
    /opt/venv/bin/pip install --upgrade pip && /opt/venv/bin/pip install jupyterlab
ENV PATH="/opt/venv/bin:$PATH"
```

Install R packages
```Dockerfile
# Install R packages
RUN R -e "install.packages(c('BiocManager', 'dplyr', 'ggplot2', 'data.table', 'future', 'cowplot'))" 
RUN R -e "BiocManager::install('tidyverse')" 
RUN R -e "BiocManager::install('Seurat')"
```

I have prepared a custom R package with all the function needed for this exam.
```Dockerfile
# Install custom R packages
RUN R -e "remotes::install_github('Maiolino-Au/PABMaiolinoPackage')"
```

Download the scripts and put them in a dedicated folder. This folder is AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```Dockerfile
# Download the R script from GitHub
RUN mkdir -p /Scripts && \
    cd /Scripts && \
    curl -O https://raw.githubusercontent.com/Maiolino-Au/ProgrAppBioinfo_Exam/main/Scripts/First.R
```

CMD
```Dockerfile
ENV SHELL=/bin/bash
CMD jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root --ServerApp.allow_origin='*' --ServerApp.token=''
```

### Build the Image
```
docker build -t maiolino_exam2025 .
```

### Run the container
To start the container use this script, either the [windods](Docker_container/script_windows_maiolino_exam2025.cmd) or [unix](Docker_container/script_unix_maiolino_exam2025.sh) version.

## Logic
This is the logic

# other
other\
\
python
```python
print("Hello, world!")
```
R
```R
print("Hello, world!")
```


