####################
#### TITLE:     Load the data frames with information about the sampled subjects: sampling without replacement in studies.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_FixRan/FixRanStudy.git/Imagen/ContK
#### First Modified: 12/05/2016
#### Notes:
#################

##
###############
### Notes
###############
##



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

# Number of studies (K) in the MA
K <- as.numeric(as.character(args))[4]

##
###############
### Loading and writing files
###############
##

# Load the study sample sizes
load(paste(LOCFRAME, '/StudySamSubj_',K,'.RData', sep=''))

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


# Load the run GT subjects
load(paste(LOCFRAME, '/GTSamSubj_',K,'.RData', sep=''))

# Select the subjects form this run
RunGTSubj <- GTSamSubj[GTSamSubj$run == RUN,'subjects']

# Write to txt file
cat(RunGTSubj,sep='\n',file=paste(wd,'/groupSubjects.txt',sep=''))
