####################
#### TITLE:     Prepare peak locations for CBMA
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 18/11/2014
#### Notes:
#################


##
###############
### Preparation
###############
##

# Location of the functions.R file: change if needed!
LocationHelperFunction <- '../functions.R'
source(LocationHelperFunction)

# Print which libraries are needed
print("Need R file: functions.R")


##########################################
##########################################
## OVERWRITE THE writeALELMax function!!!!
	# The reason is that I want to try out deleting empty studies!

writeALELMax <- function(numofsubj,indexofstudy,pathtofile,filename,pathtowrite,outputname,individFile,FirstColumn){
  ###------ HELP:
# numofsubj are the number of subjects in that particular study
# indexofstudy is the index (e.g. for the 8th study: index = 8)
# pathtofile is the absolute path to the file
# filename is how the input file is called
# pathtowrite is location where the file needs to be written
# outputname: how is output file called: WITHOUT THE .txt
# individfile: logical if procedure needs to write the individal files of each study yes/no (TRUE/FALSE)
    # Mostly not really needed (ALE uses one file with blank line between studies).
    # WARNING: IF TRUE, THEN WE NEED A FOR LOOP FOR EACH STUDY (S)!
# FirstColumn specifies which is the first x column of the file (the next two should be y and z)

if(length(grep(".txt",outputname))==1){stop("name of output file should not contain .txt!")}
  # Number of subjects
  S <- numofsubj
  # Index of study
  i <- indexofstudy
  # Open file
  dat <- try(read.table(paste(pathtofile,filename,sep=''),sep = "\t", header=TRUE, fill=TRUE,row.names=NULL),silent=TRUE)

  # Write to table: in ALE format
        #// Reference=MNI
        #// HCP, 2014
        #X Y Z

  # Individual files (for each study) needed?
  if(individFile==TRUE){
    # First create the file
    peakFile <- file(paste(pathtowrite,outputname,i,".txt",sep=""))
    writeLines(c("// Reference=MNI","// HCP,2014",paste("// Subjects=",S,sep="")), peakFile)
    close(peakFile)
    # Then append the data
    if(dim(dat)[1]==0){
      # Include empty row
      cat("\n", file = paste(pathtowrite,outputname,i,".txt",sep=""))
    } else {
      write.table(dat[,FirstColumn:(FirstColumn+2)], paste(pathtowrite,outputname,i,".txt",sep=""), sep = "\t", row.names=FALSE, col.names=FALSE, append=TRUE )
    }
  }

  # Combine them in one file with all the coordinates of the studies
	# Only write lines if there are local maxima!
	if(!dim(dat)[1]==0){
		cat("// Reference=MNI","\n",paste("// HCP,2014_",i,sep=""),"\n",paste("// Subjects=",S,sep=""),"\n",file = paste(pathtowrite, "ALEPeaks.txt",sep=""), sep="",append = TRUE)
		write.table(dat[,FirstColumn:(FirstColumn+2)], paste(pathtowrite, "ALEPeaks.txt",sep=""), sep = "\t", row.names=FALSE, col.names=FALSE, append=TRUE )
		cat("\n", file = paste(pathtowrite, "ALEPeaks.txt",sep=""), append = TRUE)
    }
}



########################################################################
########################################################################
# Take arguments from main_CBMA_K20.sh
args <- commandArgs(TRUE)

# Set working directory
wd <- as.character(args)[1]
setwd(wd)

# Index at which study we are (needed for taking the subjects from the array and for the indexofstudy in function below)
Index <- as.numeric(as.character(args)[2])
# Array of the amount of subjects in each study included in the meta-analysis
subjects <- as.numeric(as.character(args)[3])
# Array of paths were to find the files
pathtofile <- as.character(args)[4]
# How are the input files called
filename <- as.character(args)[5]
# Where to write the files
outputdir <- as.character(args)[6]
# How to name the output file
outputname <- as.character(args)[7]
# ALE format or own-meta-analysis
analysis <- as.character(args)[8]
# LOGICAL variable for ALE individual files yes/no
individFile <- as.character(args)[9]
# Where is first column with correct X coordinate
FirstColumn <- as.numeric(as.character(args)[10])



##
###############
### Write the peak locations
###############
##

# IF analysis==ALE, then use function writeALELMax
# Otherwise use function writeLocalMax (both are very similar).


## ------------------Help of function --------------------------
# numofsubj are the number of subjects in that particular study
# indexofstudy is the index (e.g. for the 8th study: index = 8)
# pathtofile is the absolute path to the file
# filename is how the input file is called
# pathtowrite is location where the file needs to be written
# outputname: how is output file called: WITHOUT THE .txt

# See function writeLocalMax
if(analysis=='ALE'){
	writeALELMax(subjects,Index,pathtofile,filename,outputdir,outputname,individFile,FirstColumn)
} else {
	writeLocalMax(subjects,Index,pathtofile,filename,outputdir,outputname,FirstColumn)
}








