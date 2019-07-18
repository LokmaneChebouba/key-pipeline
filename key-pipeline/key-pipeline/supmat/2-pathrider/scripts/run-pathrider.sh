#!/bin/bash

# Run information extraction from the graph and ICGC file
# ----------------------------------------------
# This file is part of the Supplementary Material of the submission entitled:
# Hepatocellular carcinoma computational models identify key protein-complexes associated to tumor progression
# Authors: Maxime Folschette, Vincent Legagneux, Arnaud Poret, Lokmane Chebouba, Carito Guziolowski and Nathalie Théret


#run pathrider 
#./supmat/2-pathrider/scripts/pathrider stream $1 $2 $3 #-blacklist ./supmat/0-diff-analysis/data/LIHC_primary_weakly_expressed_genes.txt
go run ./supmat/2-pathrider/scripts/pathrider/pathrider.go stream -blacklist $4 -out ./supmat/2-pathrider/data/out_pathrider.sif $1 ./supmat/2-pathrider/data/column_name.txt $2 > ./supmat/2-pathrider/pathrider_out.out

#filter the graph from a black list
#sh ./supmat/2-pathrider/scripts/remove_complexes.sh ./supmat/2-pathrider/data/out_pathrider.sif
grep -wvi -f $3 ./supmat/2-pathrider/data/out_pathrider.sif > ./supmat/2-pathrider/data/out-filtered.sif

