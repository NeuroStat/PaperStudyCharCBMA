Analyses
===========

#### Intermediate results

Raw data is processed and saved in *.rds* objects. This is done in the [IntermediateResults folder]().

#### Descriptive results
[DescriptiveResults.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/DescriptiveResults.R) contains the code to get the descriptive results (table 2).
#### ROC
[PlottingROC.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/PlottingROC.R) contains code to obtain all the ROC curves.
#### Overlap
[PlottingReliability.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/PlottingReliability.R) contains code for the measures of overlap.
#### Heatmaps
[Heatmaps.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/Heatmaps.R) contains the code to get the heatmap figures.
#### Between study variability
[BS_Var.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Analyses/BS_var.R) contains code to visualize the between study variability.

#### Configurations

The number of folds, K and the sample size in the test condition is listed in the table below.

| Folds        | K           | N in test condition  |
| ------------- |:-------------:| :-----:|
| 7 | 10 | 200 |
| 5 | 12 | 240 |
| 5 | 14 | 280 |
| 4 | 16 | 320 |
| 3 | 18 | 360 |
| 3 | 20 | 400 |
| 2 | 30 | 600 |
| 2 | 35 | 700 |


#### Results for appendix

Figures for the appendix are obtained in the [ProcesRawData.R]() file in the IntermediateResults folder.
This is because these figures are directly calculated using the raw results (and I am too lazy to change them to using the intermediate results).
