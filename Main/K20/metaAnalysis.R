####################
#### TITLE:     Perform effect size based meta-analysis
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 21/11/2014
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

# Library
require(oro.nifti)

# Take arguments from RefSamSubjK20
args <- commandArgs(TRUE)

# Location of the functions.R file: change if needed!
LocationHelperFunction <- '../functions.R'
source(LocationHelperFunction)

# Print which libraries are needed
print("Need packages: oro.nifti and R file functions.R")

# Set working directory
wd <- as.character(args)[1]
setwd(wd)

# Set seed
seed <- as.numeric(as.character(args)[2])
set.seed(seed)

# Which meta-analysis method
metamethod <- as.character(args)[3]

# Mask
locmask <- as.character(args)[4]
maskIMAGEN <- readNIfTI(locmask, verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

# Number of studies in the meta-analysis
N_S <- as.numeric(as.character(args)[5])

# Location of peak files
locStud <- as.character(args)[6]

# Name of peak files
nameStud <- as.character(args)[7]

# Directory of output (P-values)
outDir <- as.character(args)[8]

# Print the arguments: checkpoint
print(args)


##
###############
### Prepare maps for meta-analysis
###############
##

# Making progress
print(paste("STARTING ", metamethod, " Meta-Analysis",sep=""))

## Global Variables
	# Size of voxels in MM of template
	templateMM <- 3
	# R voxel size
	voxelSize <- 1
	# FWHM defined by Radua et al.
	FWHM <- 20
	# Dimension of template
	DIM <- c(53,63,46)


print(paste("Preparing maps for ", metamethod, " Effects Meta-Analysis",sep=""))
## Loading in files
allStud <- loadNewStud(N_S,locStud,nameStud)

## Model the Effect Sizes
	# Using function prepare
brain <- prepare(STUD = allStud, DIM = DIM, templateMM = templateMM, voxelSize = voxelSize, FWHM = FWHM,MASK = maskIMAGEN,origin=FSLVoxel,transformation='IMAGEN')

print(paste("Executed Preparation of ", metamethod, " Effects Meta-Analysis",sep=""))


##
###############
### Execute the meta-analalysis
###############
##


if(metamethod=="Fixed"){
MetaAnalysis <- metaAnFix(allStud, brain, DIM)
	ZstatFix <- Zstat <- MetaAnalysis[[2]]
	save(ZstatFix,file=paste(outDir,"ZstatFix",sep=""))
MetaAnalysis <- MetaAnalysis[[1]]
}
if(metamethod=="Random"){
MetaAnalysis <- metaAnMix(allStud, brain, DIM)
	ZstatRan <- Zstat <- MetaAnalysis[[4]]
	save(ZstatRan,file=paste(outDir,"ZstatRan",sep=""))

	varianceBS <- MetaAnalysis[[3]]
	save(varianceBS,file=paste(outDir,"VarBS",sep=""))
MetaAnalysis <- MetaAnalysis[[1]]
}

save(MetaAnalysis,file=paste(outDir,"MetaAnalysis",sep=""))
print(paste("Executed Weighting Studies in ", metamethod, " Effects Meta-Analysis",sep=""))



##
###############
### Calculate Null-Distribution
###############
##


# Amount of Permutations (whole brain permutations!!)
PERM <- 20
PermZStat <- c()

if(metamethod=="Fixed"){
	# MIND THE TYPE: HERE TYPE 1 FOR FIXED EFFECTS ANALYSIS!!!!!
	for(j in 1:PERM){
		print(paste((j/PERM)*100, "% done"))
		PermZStat <- cbind(PermZStat, CMakePermDis(brain = brain,allStud = allStud, type = 1, mask = maskIMAGEN))
	}
}
if(metamethod=="Random"){
	# MIND THE TYPE: HERE TYPE 2 FOR RANDOM EFFECTS ANALYSIS!!!!!
	for(j in 1:PERM){
		print(paste((j/PERM)*100, "% done"))
		PermZStat <- cbind(PermZStat, CMakePermDis(brain = brain,allStud = allStud, type = 2, mask = maskIMAGEN))
	}
}

	# Sort the distribution
	PermZStat <- sort(array(PermZStat, dim=c(prod(dim(PermZStat)))), decreasing=TRUE)


# Round to 4 digits
PermZStat <- round(PermZStat,4)

save(PermZStat,file=paste(outDir,"PermDistr",sep=""))

print(paste("Executed Null Distribution Calculation of ", metamethod, " Effects Meta-Analysis",sep=""))



##
###############
### Calculate P-Values
###############
##

# Array for P-Value, emperical values without zero and an indicator for the latter
PVal <- array(1,dim=c(prod(DIM)))
EmpZStat <- array(ZStat,dim=c(prod(DIM)))
EmpZStatExNull <- EmpZStat[which(EmpZStat>0)]
IndZStat <- which(EmpZStat>0)

# Calculate P-values using vectorized count approach (see functions.R)
PValExNull <- CPermPValCountVectorized(EmpZStatExNull, PermZStat)

# Now put these back at the correct position
	PVal[IndZStat] <- PValExNull

# Back to 3D dimension
	PVal <- array(PVal, dim=DIM)

save(PVal,file=paste(outDir,"PVal",sep=""))

print(paste("Executed P-Value Calculation of ", metamethod, " Effects Meta-Analysis",sep=""))
print(paste("ENDED ", metamethod, " Effects Meta-Analysis",sep=""))


