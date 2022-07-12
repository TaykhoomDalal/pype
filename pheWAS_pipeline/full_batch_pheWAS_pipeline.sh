#!/bin/bash

# if -abdomen is used, then run abd_full_pheWAS_pipeline.sh 
if [ $1 == "-abdomen" ]; then
    echo "Batching abd_full_pheWAS_pipeline.sh"
    sbatch abdomen_full_pheWAS_pipeline.sh 1006
    sbatch abdomen_full_pheWAS_pipeline.sh 1003
    sbatch abdomen_full_pheWAS_pipeline.sh 1019
    sbatch abdomen_full_pheWAS_pipeline.sh 17518
    sbatch abdomen_full_pheWAS_pipeline.sh 100081
    sbatch abdomen_full_pheWAS_pipeline.sh 1307
    exit 0
elif [ $1 == '-heart' ]; then
    echo "Batching heart_full_pheWAS_pipeline.sh"
    sbatch heart_full_pheWAS_pipeline.sh 1006
    sbatch heart_full_pheWAS_pipeline.sh 1003
    sbatch heart_full_pheWAS_pipeline.sh 1019
    sbatch heart_full_pheWAS_pipeline.sh 17518
    sbatch heart_full_pheWAS_pipeline.sh 100081
    sbatch heart_full_pheWAS_pipeline.sh 1307
    exit 0
elif [ $1 == '-eye' ]; then 
    echo "Batching eye_full_pheWAS_pipeline.sh"
    sbatch eye_full_pheWAS_pipeline.sh 1006
    sbatch eye_full_pheWAS_pipeline.sh 1003
    sbatch eye_full_pheWAS_pipeline.sh 1019
    sbatch eye_full_pheWAS_pipeline.sh 17518
    sbatch eye_full_pheWAS_pipeline.sh 100081
    sbatch eye_full_pheWAS_pipeline.sh 1307
    exit 0
fi