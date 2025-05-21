# Buid the docker from satijalab/seurat v5.0.0 
FROM satijalab/seurat:5.0.0

LABEL maintainer="aurelio.maiolino@edu.unito.it"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python, JupyterLab, notebook, IRkernel
RUN apt-get update && apt-get install -y \
    python3-pip python3-dev curl libzmq3-dev \
    && pip3 install --no-cache-dir jupyterlab notebook \
    && Rscript -e "install.packages('IRkernel', repos='https://cloud.r-project.org'); IRkernel::installspec(user = FALSE)"

# Install  missing R packages
RUN R -e "install.packages(c('BiocManager'))" \
    R -e "BiocManager::install(\"tidyverse\")"

# Install custom R packages
RUN R -e "remotes::install_github('Maiolino-Au/')"

# Organize
RUN mkdir Scripts 

# Download the R script from GitHub
RUN curl -o Scripts/script.R https://raw.githubusercontent.com/YOUR_USERNAME/YOUR_REPO/main/script.R

# 
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--ServerApp.root_dir=/"]
