# Modeling the COVID-19 Infection Rates by Regime Switching Unobserved Components Models

This repository contains all necessary code to replicate the results, figures and tables of the paper

**Haimerl, P. and Hartl T. 2022. Modeling the COVID-19-Infection Rates by Regime Switching Unobserved Components Models. ...**

If you use any parts of this code please cite: 

## Scripts and folders
1. **Replicate_paper.Rmd**: R-Notebook that guides you through all the necessary steps to replicate the paper
2. **Scripts**
    1. **Application.R**: Executes the model and collects the ouptut
    2. **Gridsearch.R**: Runs the 3-step grid search routine to derive the Maximum Likelihood estimates
    3. **Init.R**: Script loading all required packages and data (data provided by JH/CSSE: https://github.com/CSSEGISandData/COVID-19)
    4. **Plots.R**: Script that generates all figures
    5. **SE_from_Hessian.R**: Source code of the discontinued SEfromHessian function of the HelpersMG package (https://rdrr.io/cran/HelpersMG/src/R/SEfromHessian.R)
    6. **Tables.R**: Reproduces the tables in the paper
    7. **UC_Model.R**: Script holding all the Kim Filter as well as Kim Smoother recursions
3. **Gridsearch_UC_Model**
    1. Folder containing the results of the grid search routine as employed in the paper




