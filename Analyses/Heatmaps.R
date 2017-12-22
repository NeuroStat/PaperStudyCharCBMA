####################
#### TITLE:     Plot heatmpas to visualize activation reliability of the MAs.
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

# MNI152 template: will be used for masking in ALE
MNI152 <- readNIfTI('<<LOCATION_TO_MNI152_TEMPLATE>>/MNI152.nii')[,,]
	IDMNI152 <- MNI152 == 0

# Location to save plots for the paper
LocFileSave <- '~/Figures/'

# Results, in list according to K (number of studies in MA)
WDs <- list(
	'[1]' = "/<<LOCATION_OF_RESULTS>>/K10",
	'[2]' = "/<<LOCATION_OF_RESULTS>>/K20",
	'[3]' = "/<<LOCATION_OF_RESULTS>>/K35"
	)

# Specific functions from cowplot package.
# I gathered these from the cowplot Github page (https://github.com/wilkelab/cowplot).
# This is because the full package is not bug free on my machine as of 03/01/2017.
source('~/PaperStudyCharCBMA/Analyses/cowplot_functions.R')


##
###############
### Preparation
###############
##

# Choose your WD (in paper: 1, 6 and 8)
WD <- 8

# Setwd
setwd(WDs[[WD]])

# Libraries
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(png)
library(gridExtra)


# Extra functions
	# Going from thresholded cluster ALE map with the ALE values in significant voxels to binary maps
	# Option to provide a mask, if map is not already masked.
ThreshALECluster <- function(Thr_ALEMAP, mask = NA){
	# Already thresholded ALE maps: values larger than 0 are significant
	ID.thresholded <- Thr_ALEMAP > 0
	ALEcluster.thresholded <- Thr_ALEMAP
	ALEcluster.thresholded[ID.thresholded] <- 1
	ALEcluster.thresholded[!ID.thresholded] <- 0
	if(!is.na(mask)){
		IDmask <- mask == 0
		ALEcluster.thresholded[IDmask] <- 0
	}
	return(ALEcluster.thresholded)
}

# Threshold p-value maps
ThreshPVal <- function(PMAP, threshold, mask){
	# Put NA values, if there are some to 1
	PMAP[is.na(PMAP)] <- 1
	# Get ID for values under threshold
	ID.thresholded <- PMAP <= threshold
	# Put thresholded values to 1, the rest to 0
	PMAP.thresholded <- PMAP
		PMAP.thresholded[ID.thresholded] <- 1
		PMAP.thresholded[!ID.thresholded] <- 0
	# Safety measure if p-maps were not masked: values outside mask (i.e. zero valued in the mask) get zero
	IDmask <- mask == 0
	PMAP.thresholded[IDmask] <- 0
	return(PMAP.thresholded)
}

# Data frame with:
# -- WDs correspond to working directories (see ProcesRawData for this)
# -- number of runs/folds
# -- number of studies in the MA (K) and 
# -- number of subjects in the reference image.
DesignInfo <- data.frame(WDs = 1:8,
                         FOLDS = c(7, 5, 5, 4, 3, 3, 2, 2),
                         K = c(10, 12, 14, 16, 18, 20, 30, 35),
                         NSUBREF = c(200, 240, 280, 320, 360, 400, 600, 700))
print(DesignInfo)
NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
  unlist() %>% as.numeric()
NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
  unlist() %>% as.numeric()
NSUBREF <- DesignInfo %>% filter(WDs == WD) %>% select(NSUBREF) %>% 
  unlist() %>% as.numeric()


# Dimension of the data
DIM <- c(53,63,46)

# Vector of pooling methods
poolmeth <- c(
	'FixedEffect',
	'OLS',
	'MixedEffect'
	)
numpoolmeths <- length(poolmeth)


# Mask for thresholding the PMaps of fixed and random effects MA.
MASK <- readNIfTI(paste(WDs[[WD]],'/Run_1/GroupAnalysis/mask.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
# Mask for ALE
MASKALE <- readNIfTI('MNI152.nii', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

# Labels for group level models
ArtPOOLlabels <- c('OLS', 'Fixed Effects', 'Mixed Effects')

# Labels for meta-analyses
ArtLABMA <- c('Fixed Effects MA', 'Random Effects MA','ALE cFWE','ALE Uncorrected')

# Array of meta-analysis methods
metaMethods <- c('FixedEffUn', 'RanEffUn','ALEUn','ALECluster')
numMetaMethods <- length(metaMethods)

# number of columns in the data set
numCols <- numMetaMethods * numpoolmeths

# List of all the names for loading in the data
METHODS <- list(
	array(rep(metaMethods,each=numpoolmeths)),
	array(rep(poolmeth,numMetaMethods)))
	names(METHODS) <- c('MetaAnalysis', 'Pooling')

# The pairwise comparissons
combRuns <- c(1:NRUNS)
PAIRS <- t(combn(combRuns,2))
NPAIRS <- dim(PAIRS)[1]



##
###############
### Heatmap: NON-TRANSFORMED ALE, NON-TRANSFORMED F&R EF MA
###############
##

# NOTES:
	# Original dimension for fixed and random effects MA.
	# ALE in original dimensions.
	# And on MNI152 after which I do not need to cut the edges
	# Adding a t-statistic of a group analysis
	# Adding the heatmap of the reference images (group analyses)

##################
### Data Wrangling
##################

DIMMNI <- c(91, 109, 91)

# Start with reading in a high resolution Colin template: the Colin info and the array itself
MNIINFO <- readNIfTI('MNI152.nii', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)
	MNIPlot <- readNIfTI('MNI152.nii', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

# Now the flirted Colin
FlirtedMNI <- readNIfTI('flirtedMNI_To_Imagen.nii.gz', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

############
### PART ONE
############

# Initialize the numCols list with VOXEL x NRUNS matrices
BinaryGroupMaps <- replicate(n = numCols, expr = {array(NA, dim = c(prod(DIMMNI),NRUNS))}, simplify = FALSE)

# Loop over the pooling and meta-analysis methods
for(j in 1:numCols){
	# Loop over all the runs
	for(r in 1:NRUNS){
		if(grepl('ALE', METHODS[['MetaAnalysis']][j])){
			# Load in the data: ALE uncorrected using matlab code of Eickhoff
			if(METHODS[['MetaAnalysis']][j] == "ALEUn"){
				# Z-values and going to P-values
				Z.map <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/MetaAnalyses/ALE/',METHODS[['Pooling']][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				P.map <- 1-pnorm(Z.map, mean = 0, sd = 1)
				# Thresholding. Maps are already masked, hence NA to argument.
				Thresh.map <- ThreshPVal(P.map, 0.001, mask = NA)
				# Put this map in the list
				BinaryGroupMaps[[j]][,r] <- array(Thresh.map, dim = prod(DIMMNI))
					rm(Z.map, P.map, Thresh.map)

			}else{ # cluster Familiy wise error correction
				# ALE map: thresholded values
				filePath <- paste(WDs[[WD]],'/Run_',r,'/MetaAnalyses/ALE/',METHODS[['Pooling']][j],'/ALE/Results/',sep='')
				ALE.map <- readNIfTI(dir(path = filePath, pattern = '^OLS_cFWE05_001_', full = TRUE), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				# These are already thresholded, binarize the map here
				Thresh.map <- ThreshALECluster(Thr_ALEMAP = ALE.map)
				# Put this map in the list
				BinaryGroupMaps[[j]][,r] <- array(Thresh.map, dim = prod(DIMMNI))
					rm(ALE.map, Thresh.map)
			}
		}else{
			# Load in PVal of fixed or random effects MA
			load(paste(WDs[[WD]],'/Run_',r,'/MetaAnalyses/',METHODS[['MetaAnalysis']][j],'/',METHODS[['Pooling']][j],'/PVal',sep=''))

			# Significance testing at P < 0.005 as well as Z > 1
			# Load in the Z-values of the MA
				# We need the placeholder name for either fixed or random effects MA, this is the first 3 letters of METHODS[['MetaAnalysis']][j] (Fix or Ran)
				PLACEHOLDER <- METHODS[['MetaAnalysis']][j] %>% substr(1,3)
				load(paste(WDs[[WD]],'/Run_',r,'/MetaAnalyses/',METHODS[['MetaAnalysis']][j],'/',METHODS[['Pooling']][j],'/Zstat',PLACEHOLDER,sep=''))
			# The Z-values
			Z.map <- get(paste('Zstat',PLACEHOLDER,sep=''))
				# Make identifier for Z > 1
				GroupID.Z <- Z.map > 1
			# Threshold P-values at 0.005, while using correct mask
			P.map <- ThreshPVal(PVal, threshold = 0.005, mask = MASK)
			# Only select those with GroupID.Z = TRUE
			ThreshMap.tmp <- P.map
				ThreshMap.tmp[!GroupID.Z] <- 0
			Thresh.map <- ThreshMap.tmp
			# Put this map in the list: note that we only fill a part of the BinaryGroupMaps array here, as dimension of thresholded map of fixed and random effects < MNI dimension of ALE
			BinaryGroupMaps[[j]][c(1:prod(dim(Thresh.map))),r] <- array(Thresh.map, dim = prod(dim(Thresh.map)))
				rm(PVal, ThreshMap.tmp, GroupID.Z, PLACEHOLDER, P.map, Thresh.map)
			}
		}
	# Let us add combination of the names to the list, as a check
	method <- METHODS[['MetaAnalysis']][j]
	pooling <- METHODS[['Pooling']][j]
	NameCheck <- paste(method,':',pooling, sep = '')
		names(BinaryGroupMaps)[[j]] <- NameCheck
}

# Apply sum over the columns in each element of the list
HeatGroupMaps <- lapply(BinaryGroupMaps, function(x) apply(X = x, MARGIN = 1, FUN = sum))

# All 0 values need to get NA (for plotting purpose)
ZeroToNA <- function(x){
	IDzero <- x == 0
	x[IDzero] <- NA
	return(x)
}
	HeatGroupMaps <- lapply(HeatGroupMaps, FUN = ZeroToNA)

# Now put the images in array nVOXEL x numCols
OriginalHeatGroupMaps <- array(NA, dim = c(prod(dim(MNIPlot)), numCols))
for(j in 1:numCols){
	if(!grepl('ALE', METHODS[['MetaAnalysis']][j])){

		# Put the map inside OriginalHeatGroupMaps: note, it is not completely filled.
		OriginalHeatGroupMaps[c(1:prod(DIM)),j] <- array(HeatGroupMaps[[j]], dim = prod(DIM))

	}else{
		# We don't clip the array.
		ToWrangle <- array(HeatGroupMaps[[j]], dim = DIMMNI)
		OriginalHeatGroupMaps[,j] <- array(ToWrangle, dim = prod(dim(MNIPlot)))
	}
}

############
### PART TWO
############

# Initialize an array with nVOXEL x numCols
RefImages <- array(NA, dim = c(prod(DIM), NRUNS))

# Loop over the thresholded reference images
for(r in 1:NRUNS){
	# Load in thresholded reference image of this run
	RefImages[,r] <- array(readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/thresh_zstat1.nii', sep = ''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,], dim = prod(DIM))
}

SummedRefImages <- apply(RefImages, 1, sum, na.rm=TRUE)
	# Zero gets NA, for plotting purpose
	SummedRefImages[SummedRefImages==0] <- NA

##############
### PART THREE
##############


# Initialize arrays with nVOXEL x numCols
TstatGroups <- array(NA, dim = c(prod(DIM), NRUNS))
ZstatGroups <- array(NA, dim = c(prod(DIM), NRUNS))

# Loop over the reference t/Z statistic images
for(r in 1:NRUNS){
	# Load the t-statistic of the first run
	TstatGroups[,r] <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/stats/tstat1.nii',sep=''))[,,]

	# Load the z-statistic of the first run
	ZstatGroups[,r] <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/stats/zstat1.nii',sep=''))[,,]
}

# Calculate averages over the runs
TstatGroup <- apply(TstatGroups, 1, mean, na.rm=TRUE)
ZstatGroup <- apply(ZstatGroups, 1, mean, na.rm=TRUE)

##################
### Go to plotting
##################


# Now, we will create a large data frame based on OriginalHeatGroupMaps.
	# We need to paste all columns in one vector.
# Array with the x, y and z coordinates
XYZALE <- expand.grid(seq(1:dim(MNIPlot)[1]), seq(1:dim(MNIPlot)[2]), seq(1:dim(MNIPlot)[3]))
XYZ <- expand.grid(seq(1:DIM[1]), seq(1:DIM[2]), seq(1:DIM[3]))
PlotORHGroupMaps <- data.frame(
	'Template' = c(rep(c(matrix(FlirtedMNI, ncol = 1)), 6),												# Fixed and random MA
							rep(c(matrix(MNIPlot, ncol = 1)), 6)),														# ALE results
	'SummedVoxels' = as.factor(c(matrix(OriginalHeatGroupMaps[c(1:prod(DIM)),1:6], ncol = 1),
														 matrix(OriginalHeatGroupMaps[,c(7:12)], ncol = 1))),
	'alpha' = c(1),
	'x' = c(rep(XYZ[,1],6),rep(XYZALE[,1],6)),
	'y' = c(rep(XYZ[,2],6),rep(XYZALE[,2],6)),
	'z' = c(rep(XYZ[,3],6),rep(XYZALE[,3],6)),
	'Pooling' = c(rep(METHODS$Pooling[1:6], each = prod(DIM)),
								rep(METHODS$Pooling[c(7:12)], each = prod(dim(MNIPlot)))),
	'MA' = c(rep(METHODS$MetaAnalysis[1:6], each = prod(DIM)),
					 rep(METHODS$MetaAnalysis[c(7:12)], each = prod(dim(MNIPlot)))),
	'UsedTemplate' = c(rep(1, prod(DIM) * 6),rep(2, prod(dim(MNIPlot)) * 6))				#1: fixed and random effects MA, 2: ALE
	)
		PlotORHGroupMaps[is.na(PlotORHGroupMaps$SummedVoxels),'alpha'] <- 0
		PlotORHGroupMaps[PlotORHGroupMaps$Template == 0, 'SummedVoxels'] <- NA
		PlotORHGroupMaps$MA <- factor(PlotORHGroupMaps$MA, levels = c('ALECluster', 'ALEUn', 'FixedEffUn', 'RanEffUn'), labels = c('ALE cFWE', 'ALE Uncorrected', 'Fixed Effects MA','Random Effects MA'))
		PlotORHGroupMaps$Pooling <- factor(PlotORHGroupMaps$Pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = c('OLS', 'Fixed Effects', 'Mixed Effects'))

# Possible regions
REGIONS <- list(
	'Caudate' = 'Caudate',
	'PTS' = 'PTS'
	)
REGION <- REGIONS[['PTS']]

# Selection of Z coordinate
if(REGION == 'Caudate'){
	xyZ <- 43 							# MNI 12: interesting region = Caudate
	xyZLowRes <- 18					# Corresponds to xyZ in high res
}
if(REGION == 'PTS'){
	xyZ <- 62 							# MNI 50: interesting regions = supramarginal gyrus (posterior division) + superior parietal lobule + part of angular gyrus
	xyZLowRes <- 34					# Corresponds to xyZ in high res
}

ToPlotORHGroupMaps <- subset(PlotORHGroupMaps, (PlotORHGroupMaps$z == xyZ & PlotORHGroupMaps$UsedTemplate == 2) | (PlotORHGroupMaps$z == xyZLowRes & PlotORHGroupMaps$UsedTemplate == 1))

# Previous colour versions used the YlOrRd set of colours. Now we switched to Set1
GridSlices <- ggplot(ToPlotORHGroupMaps, aes(x = x, y = y)) + geom_tile(aes(fill = Template)) +
	scale_fill_gradient(limits = c(0, max(MNIPlot)), low = 'black', high='white', guide = FALSE) +
	geom_point(aes(x = x, y = y, colour = SummedVoxels, alpha = alpha), shape = 15, size = 0.75, na.rm = TRUE) +
	facet_wrap(Pooling ~ MA, scales = 'free') +
	scale_colour_manual(values = c(brewer.pal(name="Set1", n = NRUNS)), name = paste("Declared significant out of ",NRUNS," FOLDS: ", sep = ""), breaks =c(1:NRUNS)) +
	scale_alpha(guide = 'none') + theme_void() +
	theme(legend.position = 'top', legend.text = element_text(size=9)) +
	guides(colour = guide_legend(title.position="top", title.hjust = 0.5, nrow = 1, override.aes = list(size = 4)))

# Data frame with the average t-value of the reference images
PlotTValRef <- data.frame(
	'Template' = matrix(FlirtedMNI, ncol = 1),
	'TValue' = matrix(TstatGroup, ncol = 1),
	'alpha' = c(1),
	'x' = XYZ[,1],
	'y' = XYZ[,2],
	'z' = XYZ[,3]
	)
# Lower grey values are not plotted
PlotTValRef[PlotTValRef$Template < 30, 'alpha'] <- 0
PlotTValRef[PlotTValRef$Template == 0, 'TValue'] <- NA
# Due to transformation of the template, values at the edges are blurred, hence we put the low t-values to NA
IDZeroTVal <- !is.na(PlotTValRef$TValue) & PlotTValRef$TValue < 0.8
PlotTValRef[IDZeroTVal, 'TValue'] <- NA

	ToPlotTValRef <- PlotTValRef[PlotTValRef$z == xyZLowRes, ]
	ggplot(ToPlotTValRef, aes(x = x, y = y)) + geom_tile(aes(fill = Template)) +
		scale_fill_gradient(limits = c(0, max(MNIPlot)), low = 'black', high='white', guide = FALSE) +
		theme_void()
	ggsave(paste('FlirtedMNITemplateSliceZ_',xyZLowRes,'.png', sep=''), plot = last_plot())
	img <- readPNG(paste('FlirtedMNITemplateSliceZ_',xyZLowRes,'.png', sep=''))
	RefImageTVal <-	ggplot(ToPlotTValRef, aes(x = x, y = y)) +
		annotation_raster(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
		geom_tile(aes(fill = TValue)) +
		scale_fill_gradient(limits = c(0, ceiling(max(ToPlotTValRef$TValue, na.rm = TRUE))), low = 'yellow', high='red', na.value = 'transparent') +
		labs(fill='T-value') + theme_void() +
		theme(plot.title = element_text(hjust = 0.15, vjust = -0.2, size = 11))


# Data frame with the effect size, based on the Z statistic of the average over the reference images
PlotESRef <- data.frame(
	'Template' = matrix(FlirtedMNI, ncol = 1),
	'ES' = matrix(ZstatGroup, ncol = 1) / sqrt(NSUBREF),
	'alpha' = c(1),
	'x' = XYZ[,1],
	'y' = XYZ[,2],
	'z' = XYZ[,3]
	)
# Lower grey values are not plotted
PlotESRef[PlotESRef$Template < 30, 'alpha'] <- 0
PlotESRef[PlotESRef$Template == 0, 'ES'] <- NA
# Due to transformation of the template, values at the edges are blurred, hence we put the low ES to NA
IDZeroES <- !is.na(PlotESRef$ES) & PlotESRef$ES < 0.05
PlotESRef[IDZeroTVal, 'TValue'] <- NA

ToPlotESRef <- PlotESRef[PlotESRef$z == xyZLowRes, ]
RefImageES <-	ggplot(ToPlotESRef, aes(x = x, y = y)) +
	annotation_raster(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
	geom_tile(aes(fill = ES)) +
	scale_fill_gradient(limits = c(0, max(ToPlotESRef$ES, na.rm=TRUE)), low = 'yellow', high='red', na.value = 'transparent') +
	labs(fill='Effect Size') + theme_void() +
	theme(plot.title = element_text(hjust = 0.15, vjust = -0.2, size = 11))



# Data frame with the thresholded heatmap of the reference images
PlotRefImagesHeat <- data.frame(
	'SummedVoxels' = as.factor(matrix(SummedRefImages, ncol = 1)),
	'Template' = matrix(FlirtedMNI, ncol = 1),
	'alpha' = c(1),
	'x' = XYZ[,1],
	'y' = XYZ[,2],
	'z' = XYZ[,3]
	)
	PlotRefImagesHeat[is.na(PlotRefImagesHeat$SummedVoxels),'alpha'] <- 0
	PlotRefImagesHeat[PlotRefImagesHeat$Template == 0, 'SummedVoxels'] <- NA

ToPlotRefImagesHeat <- PlotRefImagesHeat[PlotRefImagesHeat$z == xyZLowRes, ]
	RefImageHeat <-
	ggplot(ToPlotRefImagesHeat, aes(x = x, y = y)) + geom_tile(aes(fill = Template)) +
		scale_fill_gradient(limits = c(0, max(FlirtedMNI)), low = 'black', high='white', guide = FALSE) +
		geom_point(aes(x = x, y = y, colour = SummedVoxels, alpha = alpha), shape = 15, size = 0.8, na.rm = TRUE) +
		scale_colour_manual(values = c(brewer.pal(name="Set1", n = NRUNS)), name = "", breaks =c(1:NRUNS)) +
		scale_alpha(guide = 'none') + theme_void() +
		theme(legend.position = 'right', legend.text = element_text(size=9)) +
		guides(colour = guide_legend(ncol = 1, override.aes = list(size = 4)))



# Now combine the plots using experimental development of cowplot
bottom_row <- plot_grid(RefImageHeat,  RefImageES, labels = c('B', 'C'), align = 'h', scale = 0.8, rel_widths = c(0.9,1.05))

quartz(width = 5.825243, height = 9.572815)
plot_grid(GridSlices, bottom_row, labels = c('A',''), ncol = 1, rel_heights = c(1,0.3), rel_widths = c(1,1))
ggsave(paste(getwd(), '/overlap_plots/heatmap_no_transform_',REGION,'_ES_I_',NRUNS,'.png', sep = ''), plot = last_plot())

# I used dev.size() after manually resizing the plot window to get the width and height
dev.size()
