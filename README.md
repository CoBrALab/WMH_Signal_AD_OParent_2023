# Description
 
This repository contains all analysis code and processed data to reproduce the main results for the preprint titled "Assessment of white matter hyperintensity severity using multimodal MRI in Alzheimer’s Disease" available here https://doi.org/10.1101/2023.01.20.524929.

Results in raw form and figures are also available for every main analysis.

All processed data is available in the two CSVs df_pvdeepswm.csv (for the periventricular / deep / superficial white matter parcellation) and df_lobar.csv (for the lobar parcellation). Only this data is needed to run all analyses.

Analysis code is organized such that every analysis contains a .R or .py file with the same name as the folder, and should be run in order to reproduce the results. In every folder, raw results are in the subfolder 'results_paper', and figures are in the subfolder 'visualization_paper'. If the code is run, it will create two new folders: 'results' and 'visualization' with the reproduced results and figures.

# Dependencies

- anaconda/miniconda3
- Python/3.9.7
- R/3.5.1
- pyls/0.01 (https://github.com/rmarkello/pyls)
- R modules: ggplot2, viridis, multcomp, tidyverse, broom, grid, Hmisc, corrplot
- Python packages: numpy, pandas, matplotlib, sklearn, math

