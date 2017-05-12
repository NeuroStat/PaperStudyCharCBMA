####################
#### TITLE:     Visualize between-study variability in random effects MA.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/Analyses
#### First Modified: 10/05/2016
#### Notes:
#################

##
###############
### Analyis specific directories
###############
##

# Reset workspace
rm(list=ls())

# Location to save plots for the paper
LocFileSave <- '~/Figures/'

# Results, in list according to K (number of studies in MA)
WDs <- list(
	'[1]' = "/<<LOCATION_OF_RESULTS>>/K10",
	'[2]' = "/<<LOCATION_OF_RESULTS>>/K20",
	'[3]' = "/<<LOCATION_OF_RESULTS>>/K35"
	)

##
###############
### Preparation
###############
##

# Set WD
WD <- 1

# Libraries
library(lattice)
library(grid)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(devtools)
library(fslr)
	# Print which libraries are needed
	print("Need packages: oro.nifti, fslr, lattice, ggplot2, reshape2, RColorBrewer, gridExtra and the functions.R file")

# Function to write nifti header that will be needed for correct transformation
WRITENIFTI <- function(target, source, dim){
	EndNIFTI <- nifti(img = target, dim = dim,
								datatype = 16)
		EndNIFTI@qform_code <- source@qform_code
		EndNIFTI@pixdim <- c(-1,3,3,3,1,1,1,1)
	return(EndNIFTI)
}

# Function to have a custom linear transformation (affine) through usage of FSL
CustomFlirt <- function(BrainMap, SourceMap, DIM, name, wd){
	# Create point of current wd
	currentWD <- getwd()

	# Make sure brain map is in correct dimension
	BMap <- array(BrainMap, dim = DIM)

	# We cannot transform NA values, hence we put them to 0. However, we save and return the indices of the zero elements
	ZeroID <- is.na(BMap)
	BMap[ZeroID] <- 0

	# Create a nifti image from it
	NiftiPlot <- WRITENIFTI(target = BMap, source = SourceMap, dim = DIM)

	# Write it to temp folder
	writeNIfTI(NiftiPlot, file = paste(wd, '/Transformation/', name, sep=''), gzipped = FALSE)

	# Use FSL flirt to transform to higher resolution, using the transformation matrix created above
	setwd(paste(wd, '/Transformation', sep=''))
	command <- paste('/usr/local/fsl/bin/flirt -in ', name,'.nii -ref Colin27_T1_seg_MNI_2x2x2.nii -init transformation -out flirted_',name,'.nii -omat flirted_',name,'_EPI_matrix -applyxfm', sep='')
	Sys.setenv(FSLOUTPUTTYPE="NIFTI")
	system(command)

	# Go back to previous WD
	setwd(currentWD)
}

# Number of runs/iterations.
NRUNS <- 7

# Number of studies in each meta-analysis
NSTUD <- 10

# Dimension of the data
DIM <- c(53,63,46)

# Vector of pooling methods
poolmeth <- c(
	'FixedEffect',
	'OLS',
	'MixedEffect'
	)
numpoolmeths <- length(poolmeth)


##
###############
### Tau square next to MA/GT maps
###############
##

# Variance between studies and the mask of each run
VarBS <- SWAvg <- array(NA, dim=c(prod(DIM),numpoolmeths,NRUNS))
maskR <- GTR <- tmapR <- array(NA, dim=c(DIM,NRUNS))
dim(maskR)

### Load in the data
for(r in 1:NRUNS){
	# At run:
	print(paste('At run ',r,sep=''))
  # Load in this run's GT mask
  maskR[,,,r] <- readNIfTI(paste(wd,'/Run_',r,'/GroupAnalysis/mask.nii.gz', sep=''))[,,,1]
	# Working directory for this random effects MA
	runWD <- paste(wd,'/Run_',r,'/MetaAnalyses/RanEffUn',sep='')

	# Load in the thresholded GT, then the t-values from which we will only take the significant voxels
	GTR[,,,r] <- readNIfTI(paste(wd,'/Run_',r,'/GroupAnalysis/thresh_zstat1.nii.gz', sep=''))[,,]
		idGTR <- GTR[,,,r] == 1
	tmapR.tmp <- readNIfTI(paste(wd,'/Run_',r,'/GroupAnalysis/stats/tstat1.nii.gz', sep=''))[,,]
		idzeroTmap <- tmapR.tmp == 0
		tmapR.tmp[!idGTR] <- 0
		tmapR.tmp[idzeroTmap] <- NA
	tmapR[,,,r] <- tmapR.tmp

	# Reading in the data
	for(j in 1:numpoolmeths){
		# Load in BVar
		load(paste(runWD,'/',poolmeth[j],'/VarBS',sep=''))
		# Mask the values
		varianceBS[matrix(idzeroTmap,ncol=1)] <- NA
		# Put it in the vector
		VarBS[,j,r] <- varianceBS

    # Now load in the P-values and Z-values of this random effects MA
		load(paste(runWD,'/',poolmeth[j],'/PVal',sep=''))
		load(paste(runWD,'/',poolmeth[j],'/ZstatRan',sep=''))
		# Load the weighted average values, we will only keep the significant voxels
		load(paste(runWD,'/',poolmeth[j],'/MetaAnalysis',sep=''))

		# Make identifier for Z > 1
		IDzval <- ZstatRan > 1
		# Identifier for p < 0.005
		IDpval <- PVal < .005

		# Identify masked voxels in meta-analysis
		IDzero <- MetaAnalysis == 0
		# Put the non-identified Z and p values to zero
		MetaAnalysis[!IDzval] <- 0
		MetaAnalysis[!IDpval] <- 0
		# Put zero to NA, as these are the masked voxels
		MetaAnalysis[IDzero] <- NA

		# Put in vector
		SWAvg[,j,r] <- MetaAnalysis

		# Remove object
		rm(varianceBS, MetaAnalysis, PVal, ZstatRan, IDzval, IDpval, IDzero)
	}
	# Remove runWD
	rm(runWD)
}

### Average over the runs
avVarbs <- apply(VarBS, c(1,2), mean)
avSWavg <- apply(SWAvg, c(1,2), mean)
avtmapR <- apply(tmapR, c(1,2,3), mean)

### Summary
summary(VarBS)
dim(SWAvg)
dim(VarBS)
dim(tmapR)

dim(avVarbs)
dim(avSWavg)
dim(avtmapR)


##
###############
### Plotting the variance between studies in quadrants
###############
##

# Notes:
## Now create a plot with 4 quadrants: 3 views of the brains
## with tau values and points of the significant voxels + hist.

# Colin template: high res
Colin <- readNIfTI('Colin27_T1_seg_MNI_2x2x2.nii.gz', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)
	ColinPlot <- readNIfTI('Colin27_T1_seg_MNI_2x2x2.nii.gz', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
# Transformation matrix to get to original template space of Colin (higher resolution for plotting purpose).
InvTransform <- solve(read.table('flirted_colin_EPI_matrix'))
	# Write it to .txt file, so we can use it in FSL
	WRITE <- FALSE
	if(isTRUE(WRITE)){
		write.table(InvTransform, file = paste(wd, '/Transformation/transformation', sep=''), quote = FALSE, col.names = FALSE, row.names = FALSE)
	}

# First of all, we are going to upsample our results to have higher resolution images
		# We can use the function CustomFlirt to achieve linear transformation
		# 1) Transforming the variance between study estimators and the weighted averages from the random effects MA
		TranWD <- wd
		VarBSName <- c('VarBS_FE', 'VarBS_OLS', 'VarBS_ME')
		WAvgName <- c('WAvg_FE', 'WAvg_OLS', 'WAvg_ME')
		if(isTRUE(WRITE)){
			for(t in 1:numpoolmeths){
				name.VarBS <- VarBSName[t]
				name.WAvg <- WAvgName[t]
				CustomFlirt(BrainMap = avVarbs[,t], SourceMap = Colin, DIM = DIM, name = name.VarBS, wd = TranWD)
				CustomFlirt(BrainMap = avSWavg[,t], SourceMap = Colin, DIM = DIM, name = name.WAvg, wd = TranWD)
			}
		}

		# 2) Transforming the GT t-map
		name.GT <- c('GT')
		if(isTRUE(WRITE)){
			CustomFlirt(BrainMap = avtmapR, SourceMap = Colin, DIM = DIM, name = name.GT, wd = TranWD)
		}

# Now read the images back in, while putting 0 to NA
PlotImag <- c(VarBSName, WAvgName, 'GT')
for(r in 1:length(PlotImag)){
	# Data to plot --> put values of 0 to NA (variance BS values close to 0 also to NA)
	PlotData <- readNIfTI(paste(TranWD, '/Transformation/flirted_', PlotImag[r], sep=''))[,,]
		if(grepl('VarBS', PlotImag[r])){
				PlotData[PlotData < 0.01 ] <- NA
		}else{
			PlotData[PlotData == 0] <- NA
		}
	# Assign the data to variables
	assign(paste('Plot_', PlotImag[r], sep=''),
				PlotData
			)
	# Remove data to plot
	rm(PlotData)
}

################################################################################
# Time to plot!
labels <- seq(0, 50, by  =  5)
ticks <- cbind(x = seq(0,1,length.out=length(labels)),
							y = seq(0,1,length.out=length(labels)),
							z = seq(0,1,length.out=length(labels)))

# The Z values that we will alter through plotting
xyZpoints <- c(15,35,50,66)

# The ultimate loop over the images we want to plot
for(r in 1:length(PlotImag)){
	# The image we want to plot:
	ToPlot <- get(paste('Plot_', PlotImag[r], sep=''))

	# Open new window
	# quartz(height = 2.4, width = 7.2)
	quartz(height = 2.4, width = 7.2, type='png' ,
		file = paste('/Figures/contour/',PlotImag[r],'.png',sep = ''), dpi = 100)

	# Layout properties
	layout(matrix(1:4, nrow = 1, ncol = 4, byrow = TRUE),
				 heights = c(217),
				 widths = c(170,170,170,170))
	# margin properties
	par(mar = c(0.1, 2, 4, 0.1))

	# For loop over the Z-dimension
	for(l in 1:4){
		xyZ <- xyZpoints[l]
		# The contours we want to add
		contourGT <- Plot_GT[,,xyZ]
			contourGT[!is.na(contourGT)] <- 1
			contourGT[is.na(contourGT)] <- 0
		image(ColinPlot[,,xyZ], col = grey(0:255 / 255),ann = FALSE, axes = FALSE)
		image(ToPlot[,,xyZ], add=TRUE, ann = FALSE, axes = FALSE)
		axis(2, pos = 0, at = ticks[,"z"], labels = labels)
		axis(3, pos = 1, at = ticks[,"x"], labels = labels)
		if(grepl('GT', PlotImag[r])) next
		contour(contourGT, drawlabels=FALSE, add=TRUE, col = '#3182bd')
		}
		# Remove object ToPlot
		rm(ToPlot)

		# Close connection
		dev.off()
}

