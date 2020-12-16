#!/usr/bin/env bash

##ADD YOUR SBATCH HEADERS BELOW!!!##
#SBATCH -n 1
#SBATCH --cpus-per-task=6
#SBATCH -o t6_%j.out
#SBATCH -e t6_%j.err
#SBATCH -J "PC_6"
#SBATCH -t 00:40:00
#SBATCH --constraint=highmem

./prostate > out_prostate.txt