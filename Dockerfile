FROM rocker/tidyverse:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libgdal-dev \
    libproj-dev \
    libudunits2-dev \
    libgeos-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libblas-dev \
    liblapack-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c(\
    'tidyverse', \
    'dplyr', \
    'ggplot2', \
    'jsonlite', \
    'ggridges', \
    'transport', \
    'patchwork', \
    'collapse', \
    'missForest', \
    'mice', \
    'lme4', \
    'lmerTest', \
    'broom.mixed', \
    'boot', \
    'parallel', \
    'ggseg'), \
    repos='https://cran.rstudio.com/', \
    dependencies=TRUE)"

RUN R -e "install.packages(c(\
    'lavaan'), \
    repos='https://cran.rstudio.com/', \
    dependencies=TRUE)"
