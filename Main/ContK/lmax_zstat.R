####################
#### TITLE:     Write lmax_zstat1_STD.txt file. Read from cluster_zstat.txt
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_FixRan/FixRanStudy.git/Imagen
#### First Modified: 12/02/2015
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
# Minimal cluster size
minClust <- as.numeric(as.character(args)[2])



##
###############
### Read in cluster_zstat.txt and process it
###############
##

# Data in cluster_zstat.txt file
dat <- try(read.table(file='cluster_zstat.txt', skip=1),silent=TRUE)
# If no local maxima are present: then length of dat == 1
	if(length(dat)==1){
		print("WARNING: no local maxima present. \n")
		file <- 'lmax_zstat1_STD.txt'
	cat(c("Cluster Index \t Voxels \t Value \t x \t y \t z"), file=paste(wd,'/',file,sep=''))
	} else{
		# Remove columns that we do not need
		dat <- dat[,c(1:6)]
		names(dat) <- c('Cluster Index', 'Voxels', 'Value', 'x', 'y', 'z')

		# Select only clusters with voxels > minClust and then remove column Voxels
		IDmin <- dat$Voxels >= minClust
		writeData <- dat[IDmin,-2]


		##
		###############
		### Write to file lmax_zstat1_STD.txt
		###############
		##


		write.table(writeData,file='lmax_zstat1_STD.txt', row.names=FALSE, sep='\t', quote=FALSE)
	}








