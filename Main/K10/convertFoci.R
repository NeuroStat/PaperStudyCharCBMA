####################
#### TITLE:     Convert all peak locations to MNI space (comming from matrix notation)
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 24/02/2015
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

# packages
library(AnalyzeFMRI)
library(oro.nifti)
library(devtools)


# Location of the functions.R file: change if needed!
LocationHelperFunction <- '../functions.R'
source(LocationHelperFunction)

# Print which libraries are needed
print("Need packages: AnalyzeFMRI, oro.nifti and R file: functions.R")

# Take arguments from main_CBMA_K10.sh
args <- commandArgs(TRUE)

# Set working directory
wd <- as.character(args)[1]
print(getwd())
setwd(wd)

# Name of file with local maxima to be converted
name <- as.character(args)[2]


##
###############
### Read and modify foci
###############
##

foci <- try(read.table(paste(wd,'/',name,sep=''),skip=1,sep='\t'),silent=TRUE)
	# If no foci are present: then length of foci == 1
	if(length(foci)==1){
		print("WARNING: no foci present. \n")
		print(foci)
		file <- paste('MNI_',name,sep='')
		cat(c("Cluster Index \t Value \t x \t y \t z"),file=paste(wd,'/',file,sep=''))
	} else{
			names(foci) <- c('Cluster Index', 'Value', 'x', 'y', 'z')

		tmpFoci <- data.frame(t(foci[,c(3:5)]))

		MNIFoci <- data.frame(t(apply(tmpFoci,2,SPMtoMNI,TM=TransMatrix)))[,-4]
			names(MNIFoci) <- c('xMNI', 'yMNI', 'zMNI')

		newFoci <- cbind(foci,MNIFoci)



		##
		###############
		### Write newFoci to file
		###############
		##

		file <- paste('MNI_',name,sep='')
		write.table(newFoci,paste(wd,'/',file,sep=''),quote=FALSE,row.names=FALSE,sep='\t')
		# End of if else statement
	}









