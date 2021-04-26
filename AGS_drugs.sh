#!/usr/bin/env bash
rm -rf ./sample_projects/physiboss_drugsim_gastric_AGS/
python3 ./addons/drugsims/physiboss_drugsim_meu.py -p gastric -d "AKT_inhibitor_VIII, BIRB0796, CT99021, Oxozeaenol, PKF118, PI103, PD0325901" -m both --levels "1" --cell_line AGS -cl True
find ./sample_projects/physiboss_drugsim_gastric_AGS/config/ -type f -exec sed -i "s:>1</total_concentration_levels>:>0</total_concentration_levels>:g;s:&#8236\;::g" {} \;
rename 's/\_1//g;s/\_0//g' $(find ./sample_projects/physiboss_drugsim_gastric_AGS/config/ -type f)
./reset_compile.sh
./physiboss_drugsim_gastric_AGS ./config/settings_AGS_PI103.xml

