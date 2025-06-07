# Programming Approaches for Bioinformatics Exam
Maiolino Aurelio - June 2025 \
Exam of the course Programming Approaches for Bioinformatics 

## Index
0. [Exam structure](#exam-structure)
1. [Project](#project)
2. [Docker](#docker)
3. [Logic](#logic)
4. [other](#other)



   
## Exam structure
The exam focuses primarily on programming skills; a deep understanding of the underlying results is not required. \
\
A new project will be assigned one month before each exam session (the exam description and data will be available on the moodle). Ideally, the project should be presented during the exam session immediately following its assignment. Projects may be presented within a maximum of two exam sessions. However, if the project is presented in the second session, a penalty of 2 points will be applied, reducing the maximum possible score to 28. \
\
The project can be implemented using R, Python, and/or Bash scripting. Its development is structured into three progressive stages: 

- **Minimal Requirements (Score 18–22)**
  - A high-level block diagram or procedural description must be included.
  - Develop a script that performs the required tasks.
  - Provide a README file explaining how to run the script.
- **Intermediate Level (Score 23–26)**
  - Package the script and all necessary dependencies into a Docker container. The data should remain on the local disk.
- **Advanced Level (27-30L)**
  - Transform the script into reusable functions or develop it as an R package.
  - Provide documentation of the functions through an R Markdown vignette.
  - If a R package is created, function documentation must follow the roxygen2 format.
  - Save the overall script and description in a Github.

During the oral discussion, the student must explain the decisions made throughout the project and justify the chosen approaches.\
\
It is understood that students may use Large Language Models (LLMs) for coding support. This is acceptable, provided that the student can justify the choice of tools, packages, and functions, and demonstrate an understanding of their purpose and functionality.

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
Use the dataset: “PBMCs 3k cells from a healthy donor"
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
Build the container from satijalab/seurat v5.0.0. This is an image which already contains R and Seurat, it is quicker compared to downloading them singularly and enshure the use of spcecific versions of the various libraries. \
R v4.2.0 \
Seurat v5.0.0
```Dockerfile
FROM satijalab/seurat:5.0.0
```
Install curl: it is needed to download the automatically scripts. 
```Dockerfile
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*
```
Install Python, Jupyter, and make R visible by jupiter.
```Dockerfile
RUN apt-get update && apt-get install -y \
    python3-pip python3-dev curl libzmq3-dev \
    && pip3 install --no-cache-dir jupyterlab notebook \
    && Rscript -e "install.packages('IRkernel', repos='https://cloud.r-project.org'); IRkernel::installspec(user = FALSE)"
```
satijalab/seurat already has many libraries installed, but BiocManager and tidyverse are not.
```Dockerfile
RUN R -e "install.packages(c('BiocManager'))" \
    R -e "BiocManager::install(\"tidyverse\")"
```
I have prepared a custom R package with all the function needed for this exam.
```Dockerfile
RUN R -e "remotes::install_github('Maiolino-Au/PABMaiolinoPackage')"
```
Download the scripts and put them in a dedicated folder. This folder is AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```Dockerfile
RUN mkdir -p /Scripts && \
    cd /Scripts && \
    curl -O https://raw.githubusercontent.com/Maiolino-Au/ProgrAppBioinfo_Exam/main/Scripts/First.R
```
CMD
```Dockerfile
ENV SHELL=/bin/bash
CMD jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root --ServerApp.allow_origin='*' --ServerApp.token=''
```


### Script
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


