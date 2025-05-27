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
this is the assigned project 

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
print("Wazzha!")
```
R
```R
print("Hello, world!")
print("Wazzha!")
```


