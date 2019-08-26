# Synaptopathy-Proteomics

Analysis of synaptosome TMT proteomics from __Shank2__, __Shank3__, 
__Syngap1__, and __Ube3a__ mouse models of neurodevelopmental disorders.

| Mouse Model | Disorder | REF |
| ---         | ---      | --- |
| Shank2      | ASD      | []    |
| Shank3      | ASD      |     |
| Syngap1     | ASD      |     |
| Ube3a       | ASD      |     |

## Project Organization
* Bin/ - executable scripts for the project.
* Code/ - the code that performs the projects analysis.
* Input/ - the raw data used as a starting point for the analysis.

Additional directories:
Several other directories are generated during the executation of the analysis:
* RData - directory for storing intermediate RData files.
* Figures - directory for figure output from scripts.
* Tables - directory for table output from scripts. 
* venv - Python3 virual environment for reproducibility of Python scripts.
Due to the potential for these directories to get really big, they are not 
managed by git (see `.gitignore`).

#### Analysis organization.
* [0a_Install_Dependencies.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/0a_Install_Dependencies.R)
* [0b_Functions.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/0b_Functions.R)
* [1a_TMT_Analysis_Cortex.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/1a_TMT_Analysis_Cortex.R)
* [1b_TMT_Analysis_Striatum.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/1b_TMT_Analysis_Striatum.R)
* [2_TMT_Analysis_Combined.R](https://github.com/twesleyb/Synaptopathy-Proteomics/blob/master/Code/2_TMT_Analysis_Combined.R)

## Supplemental Information.
See detailed overview of the analysis, [here](https://github.com/twesleyb/Synaptopathy-Proteomics/tree/master/Manuscript/Supplement).

## System Info:
The analysis was done using RStudio on Windows 10.

