####################
#### TITLE:     Sample subjects in K = 20 studies.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 05/09/2016
#### Notes:
#################

##
###############
### Notes
###############
##

# We have to sample 400 subjects each time.
# The subjects in the smaller groups are sampled without replacement.
# This means they cannot end up twice in the studies.

# We will sample subjects and save in data frame.
# This should be loaded in the main_CBMA_K20.sh file.


##
###############
### Preparation
###############
##

# Reset workspace
rm(list=ls())

# Set working directory
wd <- 'PaperStudyCharCBMA/Main/K20'
setwd(wd)

# Set the seed
seed <- 11121990
set.seed(seed)

# Subject ID and scanning site information:
	# OwnID is the identifier of the subject.
IDLOC <- '//<<RESTRICTED_ACCES>>/IDs'
load(IDLOC)

# Number of studies in the meta-analysis
NS <- 20

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

# Note: this is only done one time. Then the same distribution can be used in different runs with different subjects.
	# K studies with mean sample size of MeanDis and SD of SDDis
	# Have to sum to groupTotal
	# Cannot be smaller than 10 subjects per study (due to too low DOF's ==> has to be larger than 1000 for FE and ME pooling in this study)
distr <- c()
	while(sum(distr)!=groupTotal  || any(c(1:9) %in% distr){
		distr <- round(rnorm(NS,MeanDis,SDDis))
	}

if(isTRUE(WRITE)){
	write(distr,file=paste(wd,
		'/distrSampleSizes.txt',sep=''),ncolumns=NS)
} else{
	distr <- as.numeric(read.table(file=
		paste(wd,'/distrSampleSizes.txt',sep=''),header=FALSE))
}

# Sort the distribution
distr <- sort(distr)

# Number of runs that can be deployed
NRUNS <- floor(dim(OurIDs)[1]/sum(distr))


##
###############
### Sampling
###############
##

# Re-initialize the seed, to rule out effect of choosing WRITE or LOAD for
# distribution of sample sizes.
set.seed(seed)

# Initialize vectors
subjects <- study <- run <-  c()

# For loop over runs and number of studies
for(r in 1:NRUNS){
	for(i in 1:NS){
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
StudySamSubj <- data.frame('subjects' = subjects, 'study' = study, 'run' = run)

# Load the IDs again to sample the subjects in the group analyses
load(IDLOC)

# Initialize vectors
Refsubjects <- Refrun <- c()

# Sample subjects for Reference
for(r in 1:NRUNS){
	# Load the subjects from the studies, these cannot be used in the Ref
	IDstudy <- OurIDs[,'OwnID'] %in% StudySamSubj[StudySamSubj$run == r,'subjects']
	toSample <- OurIDs[!IDstudy,'OwnID']

	# Sample the data
	sampledData <- sample(toSample, groupTotal, replace = FALSE)

	# Gather in vector: subjects and run
	Refsubjects <- c(Refsubjects, sampledData)
	Refrun <- c(Refrun, rep(r, groupTotal))

	# Remove
	rm(IDstudy, toSample, sampledData)
}

# Gather in data.frame
RefSamSubj <- data.frame('subjects' = Refsubjects, 'run' = Refrun)

# If we want to check whether subjects in reference are completely separete from those in studies, then run this code
			# If FALSE, then OK
CHECK <- TRUE
if(isTRUE(CHECK)){
	for(i in 1:NRUNS){
		print(any(RefSamSubj[RefSamSubj[,'run'] == i,'subjects'] %in% StudySamSubj[StudySamSubj[,'run'] == i,'subjects']))
		print(any(StudySamSubj[StudySamSubj[,'run'] == i,'subjects'] %in% RefSamSubj[GTSamSubj[,'run'] == i,'subjects']))
	}
}

##
###############
### Save the two data frames
###############
##

# Subjects in individual studies
save(StudySamSubj, file = paste(wd,'/StudySamSubjK20.RData',sep=''))

# Subjects in the Ref
save(RefSamSubj, file = paste(wd,'/RefSamSubjK20.RData',sep=''))










