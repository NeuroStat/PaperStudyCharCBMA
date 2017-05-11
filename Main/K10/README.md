# README

## Sampling
The file [SamplingK10.R](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Main/K10/SamplingK10.R) is used to sample subjects into studies and reference images for each (unique) run/iteration.
Anonymized sampled participant IDs are saved in **R** data frames which are used later on. This file also produces the distribution of the sample sizes of all individual studies (distrSampleSizes.txt).

## Test and evaluation condition
Creating the test and evaluation condition is done with the [main_CBMA_K10.sh](https://github.com/NeuroStat/PaperStudyCharCBMA/blob/master/Main/K10/main_CBMA_K10.sh) file.
This file needs to be run for 2 times in which argument `$RUN` at line **99** should be altered. All other files are helper files and self-explanatory through the comments within the code.

