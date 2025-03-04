# Double-Spatial-Regression

Accompanying code for "Two-Stage Estimators for Spatial Confounding with Point-Referenced Data" by Nate Wiecha, Jane A. Hoppin, and Brian J. Reich.

The paper proposes "Double Spatial Regression" (DSR), aka double machine learning (Chernozhukov et al, 2018) partially linear regression applied to spatial data with Gaussian Process regression to estimate the spatial surfaces. We also extend the double machine learning (DML) partially linear model to allow vector treatment variables for greater applicability to environmental health, etc., studies. The basic idea of DSR is to subtract latent functions of space from both the treatment and response variables before estimating regression coefficients. This is helpful if both treatment and response variables exhibit spatially correlated residuals, as in that case, it is possible for standard spatial regression models to suffer from "spatial confounding," which can be thought of as an unobserved confounding variable which is spatially structured.

If you are reading this readme as part of the supplemental materials posted with the publication, consider checking the associated public GitHub repository, which has a link provided in the main paper, for any updates to the code.

The code provided allows use of DSR in applied studies, and reproduction of simulation studies provided in the paper, intended to be run on a high performance cluster. There is also code provided for the real data analysis, although data is not provided in order to protect privacy of GenX Cohort Study participants. However, we hope this code will provide a useful example for using DSR in additional real data analyses.

code//dml_function_vector.R provides code for several DSR estimators from the paper and supplement, using cross-fitting to estimate spatial trends, with a vector treatment variable. The primary function is dml_fit(), which performs DSR estimation as presented in our main paper, our recommended estimator. code//dml_function_vector_nosplit.R provides similar functions for DSR estimators without sample splitting.

code//fayetteville_analysis.R provides the code used for the real data analysis, including explanatory comments.

code//HPC includes all files needed to reproduce simulation study results. The simulation scripts in the "sim scripts" folder, as well as the submission files in the "submission scripts" folder, should be copied into the same folder, which should include an "outputs" folder, and the submission scripts should be modified as needed and submitted to the cluster. 

Within the "HPC" folder, dml_function_vector_hpc.R and dml_function_vector_nosplit_hpc.R are similar to the scripts similarly named not in the "HPC" folder, but with some functions renamed. predict_functions_hpc.R provides auxiliary functions used for either other estimators, such as the spatial linear mixed model (LMM) or gSEM for the simulation study, or functions to estimate the latent functions of space for DSR. TV_SVM_function_hpc.R provides a function for estimating the latent functions of space based on Eberts and Steinwart (2013) as mentioned in the main paper and described in detail in the supplement. dml_simulation_function_foreach_hpc.R provides a function for running a simulation given scenario parameters. dml_simulation_function_foreach_determ_hpc.R is similar but for the scenarios where the latent functions of space are deterministic, rather than drawn from Gaussian Process priors. Finally, simulation_functions_hpc.R provides functions for running simulations for the specific scenarios used in the paper's simulation study.


DoubleML (Bach et al, 2024) is the official software implementation of Chernozhukov et al (2018).

Bach P, Kurz MS, Chernozhukov V, Spindler M, Klaassen S (2024). “DoubleML: An Object-Oriented Implementation of Double Machine Learning in R.” Journal of Statistical Software, 108(3), 1–56. doi:10.18637/jss.v108.i03. 

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018). "Double/debiased machine learning for treatment and structural parameters." The Econometrics Journal, 21: C1-C68, doi:10.1111/ectj.12097.
