# Synaptopathy-Proteomics
Analysis of synaptosome TMT proteomics from __Shank2__, __Shank3__, 
__Syngap1__, and __Ube3a__ mouse models of neurodevelopmental disorders.

Main directories:
* Bin - executable scripts for the project.
* Code - the code that performs the projects analysis.
* Input - the input data.

Additional directories:
Several other directories are generated during the executation of the analysis:
* RData - directory for storing intermediate RData files.
* Figures - directory for figure output from scripts.
* Tables - directory for table output from scripts. 
* venv - Python3 virual environment for reproducibility of Python scripts.
Due to the potential for these directories to get really big, they are not 
managed by git (see `.gitignore`).

System Info:
The analysis was done using RStudio on Windows 10.

* [0a_Install-Dependencies.R]()
* [0b_Functions.R]()
* 1a_TMT_Analysis_Cortex.R
* 1b_TMT_Analysis_Striatum.R
* 2_TMT_Analysis_Combined.R

