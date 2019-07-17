#!/bin/bash

# Run information extraction from the graph and ICGC file
# ----------------------------------------------
# This file is part of the Supplementary Material of the submission entitled:
# Hepatocellular carcinoma computational models identify key protein-complexes associated to tumor progression
# Authors: Maxime Folschette, Vincent Legagneux, Arnaud Poret, Lokmane Chebouba, Carito Guziolowski and Nathalie Théret

#construction of the observation from the icgc file
sh ./supmat/1-graph-extraction/scripts/obs_construction.sh $1


