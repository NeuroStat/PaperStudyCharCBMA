####################
#### TITLE:     Create the design.mat, design.con and design.grp files for group level analyses
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 14/11/2014
#### Notes:
#################

##
###############
### Notes
###############
##

# We also need a mask, so we read in a contrast file and put all non zero elements to 1

##
###############
### Preparation
###############
##

# Library
library(oro.nifti)
	# Print which libraries are needed
	print("Need package: oro.nifti")

# Take arguments from main_CBMA_K10.sh
args <- commandArgs(TRUE)

# Set working directory
wd <- as.character(args)[1]
setwd(wd)

# Number of subjects in the study
NumSub <- as.numeric(as.character(args)[2])


##
###############
### Write design.mat
###############
##

# File connection
fileCon <- paste(wd,"design.mat",sep="")
# Text to be written to the file
cat('/NumWaves\t1
/NumPoints\t',paste(NumSub,sep=''),'
/PPheights\t\t1.000000e+00

/Matrix
',rep("1.000000e+00\n",NumSub),file=fileCon)


##
###############
### Write design.con
###############
##

fileCon <- file("design.con")
	writeLines('/ContrastName1	Group Average
/NumWaves	1
/NumContrasts	1
/PPheights		1.000000e+00
/RequiredEffect		5.034

/Matrix
1.000000e+00
',fileCon)
close(fileCon)


##
###############
### Write design.grp
###############
##

# File connection
fileCon <- paste(wd,"design.grp",sep="")
# Text to be written to the file
cat('/NumWaves\t1
/NumPoints\t',paste(NumSub,sep=''),'

/Matrix
',rep("1\n",NumSub),file=fileCon)



