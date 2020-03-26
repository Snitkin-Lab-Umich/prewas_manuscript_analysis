# prewas manuscript analysis

This repository contains the data and R scripts in support of the manuscript: [prewas: Data pre-processing for more informative bacterial GWAS](https://www.biorxiv.org/content/10.1101/2019.12.20.873158v1). The focus of this work is to assess the potential impact of variant pre-processing choices on bGWAS results.  
  
### Mansucript Authors
Katie Saund* (https://orcid.org/0000-0002-6214-6713), Zena Lapp* (https://orcid.org/0000-0003-4674-2176), Stephanie N. Thiede* (https://orcid.org/0000-0003-0173-4324), Ali Pirani (https://orcid.org/0000-0001-7810-0982), and Evan S. Snitkin (https://orcid.org/0000-0001-8409-278X)

\*Equal contribution 

### Repository contents
This repository includes the R code necessary to perform analyses and generate the figures in the manuscript.  

#### How to use this repository
##### data
The `data` direcotry contains three subdirectories: 

* `hpc`
* `local` 
* `key`

`data/hpc` contains data we generated using our high performance computing cluster (hpc). There are scripts available so that you could adapt them for your computer system, but they are slow and/or computationally intensive analyses. This directory contains some zipped files (.gz) that need to be unzipped for scripts in `lib/local` to run without error.

`data/local` contains the data the can be generated quickly on a desktop computer in R using the scripts in `lib/local` starting from data in `data/hpc` and `data/key`

Both `data/local` and `data/hpc` are subdivided by analysis. 

`data/key` contains several files necessary for plotting the data correctly (color palettes and data labels). 
##### lib
The `lib` directory contains two subdirectories: 

* `hpc`
* `local` 

`lib/hpc` contains example scripts and functions used to generate the data in `data/hpc.` These scripts will not run "as is." They are provided so that users could adapt the code to their particular computer system. 

`lib/local` contains scripts to perform any data analysis necessary to convert data in `data/hpc` into a form ready to be plotted. The script `lib/local/plot_figures.R` will use the provided data in `data/` to generate the figures in `figures/.` Scripts are written to be run from the `prewas_manuscript_analysis/` directory.

Both `lib/local` and `lib/hpc` are subdivided by analysis. 

##### figures
The plots found here were generated with the script `lib/local/plot_figures.R` The plots in this directory were finalized in Adobe Illustrator (joining panels together into one figure, resizing, etc..). 

### Sequence data  
All genome sequences available on NCBI (see Table S1). 
  
### Contributors    
Katie Saund, Zena Lapp, and Stephanie N. Thiede all contributed code to this repository. 

### Citation
prewas: Data pre-processing for more informative bacterial GWAS  
Katie Saund, Zena Lapp, Stephanie N. Thiede, Ali Pirani, Evan S Snitkin  
bioRxiv 2019.12.20.873158; doi: https://doi.org/10.1101/2019.12.20.873158  
