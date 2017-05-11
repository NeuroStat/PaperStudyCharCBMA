####################
#### TITLE:     Write Z threshold from P threshold to FDRZtresh.txt file
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 12/02/2015
#### Notes:
#################

##
###############
### Notes
###############
##


# Intermediate step needed for thresholding the group analyses.

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


##
###############
### Conversion
###############
##

Pthresh <- read.table('FDRPtresh.txt',skip=1)
revP <- 1-as.numeric(Pthresh)

Zthresh <- data.frame(qnorm(revP,0,1))
# Normally the next option is not needed as you always will have a P-value > 0
if(is.infinite(Zthresh[1,1])) Zthresh <- data.frame(100)
names(Zthresh) <- 'Z'
write.table(Zthresh, 'FDRZtresh.txt',row.names=FALSE,quote=FALSE)


# If needed to go from t to Z
# Zmap <- (sign(tmap)*(-1))*qnorm(pt(-abs(tmap),df),0,1)

