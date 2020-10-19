# PhysiBoSS 2.0: a sustainable integration of stochastic Boolean and agent-based modelling frameworks

**Version:** 2.0

**Release date:** XXX 2020

## Overview: 
PhysiBoSS 2.0 is a redesign and reimplementation of PhysiBoSS ([doi:10.1093/bioinformatics/bty766](https://doi.org/10.1093/bioinformatics/bty766)). It has been conceived as an add-on that expands the PhysiCell ([doi:10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)) agent-based functionalities with intracellular cell signalling using MaBoSS having a decoupled, maintainable and model-agnostic design. PhysiBoSS 2.0 reproduces simulations reported in the original PhysiBoSS publications and can be used with other Boolean models, for instance to predict drug synergy in a gastric adenocarcinoma cell line.

**Reference:** XXX

### How to run a PhysiBoSS sample_project inside PhysiCell:
~~~bash
git clone https://github.com/bsc-life/PhysiBoSSv2.git
cd PhysiBoSSv2
make spheroid-TNF-v2
make 
./spheroid-TNF-v2
~~~

### Key makefile rules, from [PhysiCell repository](https://github.com/MathCancer/PhysiCell):

**make \[project-name\]**: populates the indicated sample project. 
                     Use "make" to compile it. 

  PhysiBoSS \[project-name\] choices:
    drug-AGS
    spheroid-TNF
    spheroid-TNF-v2
    template-BM
    
  PhysiCell \[project-name\] choices:
    biorobots-sample
    cancer-biorobots-sample
    cancer-immune-sample 
    heterogeneity-sample
    template
    template2D 
    template3D
    virus-macrophage-sample

**make**               : compiles the current project. If no 
                     project has been defined, it first 
                     populates the cancer heterogeneity 2D 
                     sample project and compiles it 

**make list-projects** : list all available sample projects 

**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**make data-cleanup**  : clears out all simulation data 

**make reset**         : de-populates the sample project and returns to the original PhysiCell state. Use this when switching to a new PhysiCell sample project. 

Visit http://MathCancer.org/blog for tutorials and help on PhysiCell.

