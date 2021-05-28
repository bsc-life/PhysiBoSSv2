#!/bin/bash

#SBATCH --job-name="drug-prep"
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 00:05:00

module load python/3.6.1

python3 addons/drugsims/physiboss_drugsim.py -p prostate --cell_line LNCaP -d "Ipatasertib, Afatinib, Ulixertinib, Luminespib, Selumetinib, Pictilisib" -m both --levels 5 -cl True
