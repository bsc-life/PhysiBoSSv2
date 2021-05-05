#!/usr/bin/env bash
rm -rf ./sample_projects/physiboss_drugsim_gastric_AGS/
python3 ./addons/drugsims/physiboss_drugsim_meu.py -p gastric -d "AKT_inhibitor_VIII, BIRB0796, CT99021, Oxozeaenol, PKF118, PI103, PD0325901" -m both --levels "0" --cell_line AGS -cl True -s "./sample_projects/gastric/config/gastric_drug_sensitivity2.csv"
for f in ./sample_projects/physiboss_drugsim_gastric_AGS/config/*.xml; do gawk -i inplace '/multiplier_reporter/ && f++>0 {next} 1' $f; done
./reset_compile.sh
./physiboss_drugsim_gastric_AGS ./config/settings_AGS_PI103.xml
