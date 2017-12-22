####################
#### TITLE:     Sample subjects in K studies to check effect of adding more studies: continuously increase K.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_FixRan/FixRanStudy.git/Imagen/35Studies
#### First Modified: 20/11/2017
#### Notes:
#################

##
###############
### Notes
###############
##

# The subjects in the smaller groups are sampled without replacement.

# We will sample subjects and save in data frame.
# This should be loaded in the master file.


##
###############
### Preparation
###############
##

# Reset workspace
rm(list=ls())

# Set working directory
wd <- '/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_FixRan/FixRanStudyGit.git/Imagen/ContK'
setwd(wd)

# Set the seed
seed <- 11121990
set.seed(seed)

# Scanning site information
IDLOC <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/OurIDs'
load(IDLOC)

# As information, the next subjects are deleted from the sampling process:
			# We do not remove sujbects here anymore, this is already done in the OurIDs file!
			# I just let them here as a reference!!!!
		remove <- c(1284,680,814,1460,1470,788,																															# Corrupted files, visual inspection
		40, 54, 61, 69, 155, 349, 384, 591, 869, 1219, 																											# Less amount of total masked voxels
		43,77,89,169,180,235,375,394,471,489,544,637,744,787,903,957,1043,1065,1129,1168,1232,1324,1335,1351,																					# Less amount of total masked voxels
		10, 166, 245, 373, 428, 474, 482, 491, 536, 565, 588, 600, 613, 618, 662, 668, 685, 699, 722, 747, 799, 829, 890, 987, 997, 1013, 1078, 1097,	# Less amount of total masked voxels
		1140, 1212, 1217, 1228, 1302, 1310, 1368, 1403, 1406, 1462, 1478,																		# No activity in upper brain regions
		68, 250, 328, 475, 551, 757, 1085, 1420)																														# No activity in upper brain regions

# Number of studies in the meta-analysis
NS <- 10:35

# Mean and standard deviation of distribution of sample sizes
MeanDis <- 20
SDDis <- 5

# Total amount of subjects over all studies
groupTotal <- MeanDis * NS

##
###############
### Creating the distribution of the sample sizes in the individual studies
###############
##

# If this file is run for the first time, then write the distribution of the
	# sample sizes. Otherwise load it.
WRITE <- FALSE

# Run over the amount of studies in the MA
for(s in 1:length(NS)){
    # Note: this is only done once for each setting. Then the same distribution can be used in different runs with different subjects.
  	# K studies with mean sample size of MeanDis and SD of SDDis
  	# Have to sum to groupTotal[s]
  	# Cannot be smaller than 10 subjects per study (due to too low DOF's ==> has to be larger than 1000 for FE and ME pooling in this study)
  distr <- c()
  	while(sum(distr)!=groupTotal[s] || any(c(1:9) %in% distr)){
  		distr <- round(rnorm(NS[s],MeanDis,SDDis))
  	}
  
  if(isTRUE(WRITE)){
  	write(distr,file=paste(wd,
  		'/_K_',9 + s,'/distrSampleSizes_',9 + s,'.txt',sep=''),ncolumns=NS[s])
  } else{
  	distr <- as.numeric(read.table(file=
  		paste(wd,'/_K_',9 + s,'/distrSampleSizes_',9 + s,'.txt', sep = ''),header=FALSE))
  }
}
  


##
###############
### Sampling
###############
##

# Re-initialize the seed, to rule out effect of choosing WRITE or LOAD for
# distribution of sample sizes.
set.seed(seed)

# Assemble the number of runs 
NRUNS_vector <- NTEST <- c()

# For loop over the amount of studies in the MA
for(s in 1:length(NS)){
  # Read in the distribution
  distr <- as.numeric(read.table(file=
       paste(wd,'/_K_',9 + s,'/distrSampleSizes_',9 + s,'.txt', sep = ''),header=FALSE))
  
  # Sort the distribution
  distr <- sort(distr)
  
  # Number of runs that can be deployed
  NRUNS <- floor(dim(OurIDs)[1]/sum(distr))
  NRUNS_vector <- c(NRUNS_vector, NRUNS)
  
  # Sample size of test condition
  NTEST <- c(NTEST, groupTotal[s])

  
  # Initialize vectors
  subjects <- study <- run <-  c()
  
  # For loop over runs and number of studies
  for(r in 1:NRUNS){
  	for(i in 1:NS[s]){
  		toSample <- OurIDs[,'OwnID']
  
  		# Sample the data
  		sampledData <- sample(toSample,distr[i],replace=FALSE)
  
  		# Gather subjects, study and run
  		subjects <- c(subjects, sampledData)
  		study <- c(study, rep(i,distr[i]))
  		run <- c(run, rep(r, distr[i]))
  
  		# Remove that sampled data from total pool and OurIDs
  		ID.remove <- OurIDs[,'OwnID'] %in% sampledData
  		OurIDs <- OurIDs[!ID.remove,]
  
  		# Remove objects
  		rm(toSample, sampledData, ID.remove)
  	}
  }
  
  # Gather vectors in data.frame
  StudySamSubj <- data.frame('subjects' = subjects, 'study' = study, 'run' = run, 'K' = NS[s])
  
  # Load the IDs again to sample the subjects in the group analyses
  load(IDLOC)
  
  # Initialize vectors
  GTsubjects <- GTrun <- c()
  
  # Sample subjects for GT
  for(r in 1:NRUNS){
  	# Load the subjects from the studies, these cannot be used in the GT
  	IDstudy <- OurIDs[,'OwnID'] %in% StudySamSubj[StudySamSubj$run == r,'subjects']
  	toSample <- OurIDs[!IDstudy,'OwnID']
  
  	# Sample the data
  	sampledData <- sample(toSample, groupTotal[s], replace = FALSE)
  
  	# Gather in vector: subjects and run
  	GTsubjects <- c(GTsubjects, sampledData)
  	GTrun <- c(GTrun, rep(r, groupTotal[s]))
  
  	# Remove
  	rm(IDstudy, toSample, sampledData)
  }
  
  # Gather in data.frame
  GTSamSubj <- data.frame('subjects' = GTsubjects, 'run' = GTrun, 'K' = NS[s], 'GroupTotal' = groupTotal[s])
  
  # If we want to check whether subjects in GT are completely separete from those in studies, then run this code
  			# If FALSE, then OK
  CHECK <- TRUE
  if(isTRUE(CHECK)){
  	for(i in 1:NRUNS){
  		print(any(GTSamSubj[GTSamSubj[,'run'] == i,'subjects'] %in% StudySamSubj[StudySamSubj[,'run'] == i,'subjects']))
  		print(any(StudySamSubj[StudySamSubj[,'run'] == i,'subjects'] %in% GTSamSubj[GTSamSubj[,'run'] == i,'subjects']))
  	}
  }
  
  
  
  ##
  ###############
  ### Save the two data frames
  ###############
  ##
  
  # Subjects in individual studies
  #save(StudySamSubj, file = paste(wd,'/_K_', 9 + s, '/StudySamSubj_', 9 + s, '.RData',sep=''))
  
  # Subjects in the GT
  #save(GTSamSubj, file = paste(wd,'/_K_', 9 + s, '/GTSamSubj_', 9 + s, '.RData',sep=''))
  
}
  

# Print number of runs (folds), next to number of studies (K)
infoRUN_K <- data.frame('NRUNS' = NRUNS_vector, 'K' = NS, 'NTEST' = NTEST)
print(infoRUN_K)


# Select the number of studies and show the number of associated runs
KSEL <- c(10,12,14,16,18,20,30,35)
library(dplyr)
infoRUN_K %>% filter(K %in% KSEL)

# Some checks
infoRUN_K %>% mutate(TOTALN = NRUNS * K * 20,
              CHECKN = K * 20,
              CHECKN_L = ifelse(CHECKN == NTEST, TRUE, FALSE))



