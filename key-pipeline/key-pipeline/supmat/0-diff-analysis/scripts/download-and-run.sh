#!/bin/bash

# Fetch and filter data, and perform differential expression analysis and clustering
# ----------------------------------------------------------------------------------
# This file is part of the Supplementary Material of the submission entitled:
# Hepatocellular carcinoma computational models identify key protein-complexes associated to tumor progression
# Authors: Maxime Folschette, Vincent Legagneux, Arnaud Poret, Lokmane Chebouba, Carito Guziolowski and Nathalie Théret

# This script downloads ICGC files for differential expression analysis and calls
# data filtering (dataset_filtering.sh) and differential expression analysis
# (diffexp_and_clustering.R) scripts.
#
# Install R and required packages with:
#   $ conda install r r-gplots r-RColorBrewer
#
# Warning:
#   The donwloads may take some time, but they are only performed once.
#
# Usage:
#   bash download-and-run.sh
#

if [ ! -f sample.LIHC-US.tsv ] && [ ! -f specimen.LIHC-US.tsv ] && [ ! -f exp_seq.LIHC-US.tsv ]
then
  wget https://dcc.icgc.org/api/v1/download?fn=/release_21/Projects/LIHC-US/sample.LIHC-US.tsv.gz -O sample.LIHC-US.tsv.gz 
  wget https://dcc.icgc.org/api/v1/download?fn=/release_21/Projects/LIHC-US/specimen.LIHC-US.tsv.gz -O specimen.LIHC-US.tsv.gz
  wget https://dcc.icgc.org/api/v1/download?fn=/release_21/Projects/LIHC-US/exp_seq.LIHC-US.tsv.gz -O exp_seq.LIHC-US.tsv.gz
  gunzip -k sample.LIHC-US.tsv.gz specimen.LIHC-US.tsv.gz exp_seq.LIHC-US.tsv.gz
fi

bash dataset_filtering.sh
R -f diffexp_and_clustering.R

