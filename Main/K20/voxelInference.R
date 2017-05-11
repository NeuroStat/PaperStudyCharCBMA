####################
#### TITLE:     Perform voxelwise inference for group analysis
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 10/02/2015
#### Notes: FDR correction, Benjamin & Hochberg procedure
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
library(AnalyzeFMRI)
library(lattice)
library(oro.nifti)
	# Print which libraries are needed
	print("Need packages: AnalyzeFMRI, lattice and oro.nifti")

# Take arguments from main_CBMA_K20
args <- commandArgs(TRUE)

# Set working directory
wd <- as.character(args)[1]
setwd(wd)
# Directory of pval
pmapDir <- as.character(args)[2]
pmap <- readNIfTI(paste(pmapDir,"/pval.nii.gz",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

# Mask: to get proper amount of null hypotheses
mask <- readNIfTI(paste(wd,"/mask.nii.gz",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]

# Significance level
signLevel <- as.numeric(as.character(args)[3])


##
###############
### FDR Correction
###############
##

# Dimension of the pmap
DIM <- dim(pmap)

## Mask the p-values map
idMask <- mask==0
MPmap <- pmap
MPmap[idMask] <- NA

## Adjust P-values with Benjamin & Hochberg procedure (built in R function)
adjPvals <- array(p.adjust(MPmap,method="BH"),dim=DIM)

## Threshold at signLevel BH p-value: significant P-values get 1!
idP <- adjPvals <= signLevel
threshPval <- adjPvals
threshPval[!idP] <- 0
threshPval[idP] <- 1


##
###############
### Make a levelplot
###############
##

switchAxis <- threshPval[c(dim(threshPval)[1]:1),,]
pdf('GroupAnalysis.pdf')
levelplot(threshPval,main=paste('Group Analysis: BH adjusted p-values < ', signLevel,sep=''),
	pretty=TRUE, col.regions = terrain.colors, xlab = 'X', ylab = 'Y',
	xlim=c(0,DIM[1]),ylim=c(0,DIM[2]))
dev.off()



##
###############
### Make nifti file
###############
##
idNA <- is.na(threshPval)
threshPval[idNA] <- 0

niftiimage <- nifti(img=threshPval,dim=DIM)
writeNIfTI(niftiimage,filename=paste(wd,'/thresh_zstat1',sep=''),gzipped=FALSE)











