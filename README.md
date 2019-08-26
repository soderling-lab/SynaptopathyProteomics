# Synaptopathy-Proteomics

Analysis of synaptosome TMT proteomics from __Shank2__, __Shank3__, 
__Syngap1__, and __Ube3a__ mouse models of neurodevelopmental disorders.

| Mouse Model | Disorder | Reference |
| ---         | ---      | --- |
| Shank2      | ASD      | [REF](url) |
| Shank3      | ASD      | [Wang et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27161151) |
| Syngap1     | ID, ASD, epilepsy |[Kim et al., 2003](https://www.ncbi.nlm.nih.gov/pubmed/12598599) |
| Ube3a       | Angelman Syndrome | [Jiang et al. 1998,](https://www.ncbi.nlm.nih.gov/pubmed/9808466) |

## Project Organization
* bin/ - executable scripts for the project.
* code/ - the code that performs the projects analysis.
* functions/ - functions used in the anallysis.
* input/ - the raw data used as a starting point for the analysis.
* tables/ - directory for table output from scripts. 

## Additional directories:
Several other directories are used to store oup:
* data - directory for storing intermediate data files.
* figures - directory for figure output from scripts.
Due to the potential for these directories to get really big, they are not 
managed by git (see `.gitignore`).

## Analysis organization.
* [0a_Install_Dependencies.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/0a_Install_Dependencies.R)
* [0b_Functions.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/0b_Functions.R)
* [1a_TMT_Analysis_Cortex.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/1a_TMT_Analysis_Cortex.R)
* [1b_TMT_Analysis_Striatum.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/1b_TMT_Analysis_Striatum.R)
* [2_TMT_Analysis_Combined.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/2_TMT_Analysis_Combined.R)

## Supplemental Information.
For a detailed overview of tha analysis, see [here](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/Manuscript/Supplement).

#### System Info:
The analysis was done using the Windows Substem for Linux (WSL) on Windows 10.  
sysname        "Windows"  
release        "10 x64"          
version        "build 18362"     
nodename       "WINDOWS10"       
machine        "x86-64"          
login          "Tyler Bradshaw"  
user           "Tyler Bradshaw"  
effective_user "Tyler Bradshaw"  
