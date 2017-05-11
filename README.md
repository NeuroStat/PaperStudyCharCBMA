The influence of study characteristics on coordinate-based fMRI meta-analyses.
===========

## Scripts
For each set size (K), there is a folder containing a file called *master_fixran_HPC.sh*. This file contains the main body of the design of the study (figure 2 in the paper).
This is run on the UGent HPC infrastructure.  
We obtain with this:
  - Group analysis of evaluation condition (eg. N=200)
  - Analyses of small studies using 3 pooling methods at second phase GLM
  - Fixed and random effects meta-analyses using *R* based on results from above.

The analyses for Activation Likelihood Estimation were performed using MATLAB scripts from prof. dr. Simon Eickhoff. 


## Analyses
PlottingROC.R contains code to obtain all the ROC curves. PlottingReliability.R contains code for the descriptive results, the measures of overlap and the heatmaps.
PlottingBS_Var.R contains code to visualise the between study variability.