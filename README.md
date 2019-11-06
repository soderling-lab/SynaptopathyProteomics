# Synaptopathy-Proteomics

Analysis of synaptosome TMT proteomics from __Shank2__, __Shank3__, 
__Syngap1__, and __Ube3a__ mouse models of human autism spectrum disorders 
(ASD), intellectual disability (ID), and epilepsy. 

-------------------------------------------------------------------------------

| Mouse Model | Disorder | Reference |
| ---         | ---      | --- |
| Shank2      | ASD      | [REF](url) |
| Shank3      | ASD      | [Wang et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27161151) |
| Syngap1     | ID, ASD, epilepsy |[Kim et al., 2003](https://www.ncbi.nlm.nih.gov/pubmed/12598599) |
| Ube3a       | Angelman Syndrome | [Jiang et al. 1998,](https://www.ncbi.nlm.nih.gov/pubmed/9808466) |

-------------------------------------------------------------------------------

## Repository Organization
### Directories
* [bin/](https://github.com/twesleyb/SynaptopathyProteomics/tree/master/bin) - executable scripts for the project.
* [analysis/](https://github.com/twesleyb/SynaptopathyProteomics/tree/master/analysis) - this directory contains scripts used to perform the data analaysis. The scripts are named in the order they should be executed in order to complete the analysis. 
* [R/](https://github.com/twesleyb/SynaptopathyProteomics/tree/master/R) - functions used in the analysis.
* [rdata/](https://github.com/twesleyb/SynaptopathyProteomics/tree/master/rdata) - directory for storing intermediate data files.
* [input/](https://github.com/twesleyb/SynaptopathyProteomics/tree/master/input) - the raw data used as a starting point for the analysis.
* [tables/](https://github.com/twesleyb/SynaptopathyProteomics/tree/master/tables) - directory for table output from scripts. 
* [data/](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/data) - directory for storing data for R package.
* [figures/](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/figures) - directory for plots and figures.

## Data Analysis 
For a detailed overview of the analysis, please see the supplemental information [here](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/code/README.md).

## Reproducibility

You can reproduce the developing environment in which the analysis was performed using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).
Use `conda env create -f venv.yml` to  to set up all R and Python dependencies.

This will create a virtual environment. 
Use `conda activate SynaptopathyProteomics` to activate it.

This work is indebted to numerous others who have shared their ideas, software,
and time. In particular, the following open source packages were essential in 
completion of this work:
* __WGCNA__: [package](https://cran.r-project.org/web/packages/WGCNA/index.html) | [publication](https://www.ncbi.nlm.nih.gov/pubmed/19114008)
* __anRichment__: [package]() 
* __EdgeR__: [package](https://bioconductor.org/packages/release/bioc/html/edgeR.html) | [publication](https://www.ncbi.nlm.nih.gov/pubmed/19910308)
* __NetRep__: [package](https://cran.rstudio.com/web/packages/NetRep/index.html) | [publication](https://www.ncbi.nlm.nih.gov/pubmed/27467248)
* __leidenalg__: [package](https://pypi.org/project/leidenalg/) | [pre-print publication](https://arxiv.org/abs/1810.08473)  

For more information about these packages, please see their respective publications. 
See also, a complete list of [R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/bin/r_requirements.txt) and 
[Python](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/bin/python_requirements.txt) dependencies. 

## System Info
The analysis was done using the Windows Substem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)) on a Windows 10 PC. 

|          |             |
| ---      | ---         |
| sysname  | Windows     |
| release  | 10 x64      |
| version  | build 18362 |
| nodename | WINDOWS10   | 
| machine  | x86-64      | 
