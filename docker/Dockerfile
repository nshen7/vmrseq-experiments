## Base image https://hub.docker.com/u/bioconductor/bioconductor_docker:RELEASE_3_15/
FROM bioconductor/bioconductor_docker:RELEASE_3_15

MAINTAINER Ning Shen

LABEL version="0.0.0"

## install CRAN packages

RUN R -e "devtools::install_version('pheatmap', version = '1.0.12', dependencies = T)"
RUN R -e "devtools::install_version('tidyverse', version = '1.3.2', dependencies = T)"
RUN R -e "devtools::install_version('data.table', version = '1.14.2', dependencies = T)"
RUN R -e "devtools::install_version('stringr', version = '1.4.0', dependencies = T)"
RUN R -e "devtools::install_version('readr', version = '2.1.2', dependencies = T)"
RUN R -e "devtools::install_version('readxl', version = '1.4.1', dependencies = T)"
RUN R -e "devtools::install_version('gamlss', version = '5.4-12', dependencies = T)"
RUN R -e "devtools::install_version('gamlss.dist', version = '6.0-5', dependencies = T)"
RUN R -e "devtools::install_version('recommenderlab', version = '1.0.3', dependencies = T)"
RUN R -e "devtools::install_version('Matrix', version = '1.5-3', dependencies = T)"
RUN R -e "devtools::install_version('colorspace', version = '2.1-0', dependencies = T)"
RUN R -e "devtools::install_version('colourvalues', version = '0.3.7', dependencies = T)"
RUN R -e "devtools::install_version('colorspace', version = '2.1-0', dependencies = T)"
RUN R -e "devtools::install_version('colorspace', version = '2.1-0', dependencies = T)"

## install Bioconductor packages

RUN R -e "BiocManager::install('biomaRt', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('limma', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('DESeq2', version='3.15', update = FALSE)"
# RUN R -e "BiocManager::install('dmrseq', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('GEOquery', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('BiocParallel', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('clusterProfiler', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('GenomicFeatures', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('annotatr', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('AnnotationDbi', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('SummarizedExperiment', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('GenomicRanges', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('DelayedMatrixStats', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('HDF5Array', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('DelayedArray', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('Biobase', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('org.Mm.eg.db', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('org.Hs.eg.db', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene', version='3.15', update = FALSE)"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19', version='3.15', update = FALSE)"

## temporary

# RUN apt-get install -y --no-install-recommends liblzma-dev

RUN R -e "BiocManager::install('MEDIPS', version='3.15', update = FALSE)"

## alias to launch rstudio server
# RUN echo "alias rstudio='/init'" >> ~/.bashrc

# ENV BASEDIR=/home/rstudio
# WORKDIR /home/rstudio

CMD ["/bin/bash"]