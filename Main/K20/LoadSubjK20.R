####################
#### TITLE:     Load the data frames with information about the sampled subjects: sampling without replacement in studies.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 12/05/2016
#### Notes:
#################

##
###############
### Notes
###############
##

# NOTE: ALL SUBJECT IDS ARE ANONYMIZED!!

# Sampling subjects has been done using the SamplingK20.R file.
# We use this file to write text files with the subject numbers
# to the StudySamSubjDouble.sh file.


# SETTING:
# The subjects in the smaller groups are sampled without replacement.
# This means they cannot end up twice in the studies.
# We will use this setting to calculate reliability of a meta-analysis.



##
###############
### Preparation
###############
##

# Take arguments from master file
args <- commandArgs(TRUE)

# Set working directory
wd <- as.character(args)[1]
setwd(wd)

# Which run are we in?
RUN <- as.numeric(as.character(args)[2])

# Location of the data frames
LOCFRAME <- as.character(args)[3]


##
###############
### Loading and writing files
###############
##

# Load the study sample sizes
load(paste(LOCFRAME, '/StudySamSubjK20.RData', sep=''))

# Select the subjects from this run
RunSubj <- StudySamSubj[StudySamSubj$run == RUN,]

# Number of studies
NS <- length(unique(StudySamSubj$study))

# For loop over the studies
for(s in 1:NS){
  # Take the subjects
  sampledData <- RunSubj[RunSubj$study == s,'subjects']
  # Write them to txt file
  cat(sampledData,file=paste(wd,"/Study_",s,"/study_",s,".txt",sep=""),sep='\n')
}


###############################
###############################


# Load the run reference subjects
load(paste(LOCFRAME, '/RefSamSubjK20.RData', sep=''))

# Select the subjects form this run
RunRefSubj <- RefSamSubj[RefSamSubj$run == RUN,'subjects']

# Write to txt file
cat(RunRefSubj,sep='\n',file=paste(wd,'/groupSubjects.txt',sep=''))












