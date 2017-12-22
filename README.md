The Influence of Study-Level Inference Models and Study Set Size on Coordinate-Based fMRI Meta-Analyses.
===========

## Scripts
Located in [Main folder](https://github.com/NeuroStat/PaperStudyCharCBMA/tree/master/Main). For each set size (K), there is a folder containing a file called *master_fixran_HPC.sh*. This file contains the main body of the design of the study (figure 2 in the paper).
This is run on the UGent HPC infrastructure.  
We obtain with this:
  - Group analyses of evaluation condition.
  - Analyses of small studies using the 3 group level models (FE, ME and OLS) at second phase GLM
  - Fixed and random effects coordinate-based meta-analyses using *R* based on results from individual studies.

The analyses for Activation Likelihood Estimation were performed using MATLAB scripts obtained from prof. dr. Simon Eickhoff.


## Analyses
Located in [Analyses folder](https://github.com/NeuroStat/PaperStudyCharCBMA/tree/master/Analyses).
Within this folder, we have a directory containing all figures from the paper. And a directory containing the intermediate results in **R** objects. 
#### Descriptive results
[DescriptiveResults.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/DescriptiveResults.R) contains the code to get the descriptive results (table 3).
#### ROC
[PlottingROC.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/PlottingROC.R) contains code to obtain all the ROC curves.
#### Overlap
[PlottingReliability.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/PlottingReliability.R) contains code for the measures of overlap.
#### Heatmaps
[Heatmaps.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/Heatmaps.R) contains the code to get the heatmap figures.
#### Between study variability
[BS_Var.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/BS_var.R) contains code to visualize the between study variability.
