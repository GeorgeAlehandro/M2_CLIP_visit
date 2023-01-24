FROM rocker/rstudio:4.2.2
 
MAINTAINER Geroge Alehandro Saad <george-alehandro.saad@protonmail.com>
 
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y \
    libhdf5-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libxt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk40 \
  && apt-get clean all && \
  apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
 
RUN Rscript -e "install.packages(c('devtools', 'phytools'));"

RUN Rscript -e "devtools::install_github("PeterZZQ/VeloSim");"
