# Synaptopathy-Proteomics

Analysis of synaptosome TMT proteomics from __Shank2__, __Shank3__, 
__Syngap1__, and __Ube3a__ mouse models of human autism spectrum disorders 
(ASD), intellectual disability (ID), and epilepsy. 

| Mouse Model | Disorder | Reference |
| ---         | ---      | --- |
| Shank2      | ASD      | [REF](url) |
| Shank3      | ASD      | [Wang et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27161151) |
| Syngap1     | ID, ASD, epilepsy |[Kim et al., 2003](https://www.ncbi.nlm.nih.gov/pubmed/12598599) |
| Ube3a       | Angelman Syndrome | [Jiang et al. 1998,](https://www.ncbi.nlm.nih.gov/pubmed/9808466) |

## Project Organization
* [bin/](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/bin) - executable scripts for the project.
* [code/](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/code) - the code that performs the projects analysis.
* [functions](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/functions) - functions used in the anallysis.
* [input/](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/input) - the raw data used as a starting point for the analysis.
* [tables/](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/input) - directory for table output from scripts. 

## Additional directories
Several other directories are used to store output generated during the analysis:
* data/ - directory for storing intermediate data files.
* figures/ - directory for figure output from scripts.
Due to the potential for these directories to get really big, they are not 
managed by git (see `.gitignore`).

## Supplemental Information
For a detailed overview of the analysis, see the supplemental information [here](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/code/README.md).

## Acknowledgements
This work is indebted to numerous others who have shared their ideas, software,
and time. In particular, the following open source packages were essential in 
completion of this work:
[WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html)
[EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) 
[NetRep](https://cran.rstudio.com/web/packages/NetRep/index.html)
[leidenalg](https://pypi.org/project/leidenalg/)

## System Info
The analysis was done using the Windows Substem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)) on Windows 10. 

|     |     |
| --- | --- |
| sysname | Windows |
| release | 10 x64 |
| version | build 18362 |
| nodename | WINDOWS10 | 
| machine | x86-64 | 
