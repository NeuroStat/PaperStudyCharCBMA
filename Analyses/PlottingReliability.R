####################
#### TITLE:     Calculate activation reliability of the MAs.
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

# Choose your WD
WD <- 3

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

# Calculate overlap between binary (thresholded) maps
OVERLAP <- function(mapA, mapB, mask, precision = 4){
	# Transform to one vector per map, if not same length, provide error
	mapA <- matrix(mapA, ncol = 1) ; mapB <- matrix(mapB, ncol = 1); mask <- matrix(mask, ncol = 1)
		if(!identical(length(mapA), length(mapB), length(mask))) stop('Map A, B and mask should be of same length!')
	# First put voxels outside mask to NA
	IDmask <- mask == 0
	mapA[IDmask] <- NA; mapB[IDmask] <- NA
	# Calculate intersection, Va and Vb
	both <- mapA + mapB
	Vab <- length(which(both == 2))
	Va <- length(which(mapA == 1))
	Vb <- length(which(mapB == 1))
	# Calculate value
	if(sum(Vab,Va,Vb) == 0){
		value <- 0
	}else{
		value <- Vab / (Va + Vb - Vab)
	}
		return(round(value, precision))
}


# Number of runs.
if(WD == 1) NRUNS <- 7
if(WD == 2) NRUNS <- 3
if(WD == 3) NRUNS <- 2

# Number of studies in each meta-analysis
if(WD == 1) NSTUD <- 10
if(WD == 2) NSTUD <- 20
if(WD == 3) NSTUD <- 35

# Number of subjects in the reference image
if(WD == 1) NSUBREF <- 400
if(WD == 2) NSUBREF <- 400
if(WD == 3) NSUBREF <- 700

# Dimension of the data
DIM <- c(53,63,46)

# Vector of pooling methods
poolmeth <- c(
	'FixedEffect',
	'OLS',
	'MixedEffect'
	)
numpoolmeths <- length(poolmeth)


# Mask for thresholding the PMaps of fixed and random effects MA: note that we used one universal mask throughout the entire study.
MASK <- readNIfTI(paste(WDs[[WD]],'/Run_1/GroupAnalysis/mask.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
# Mask for ALE
MASKALE <- readNIfTI('/Users/hanbossier/Dropbox/PhD/PhDWork/Meta Analysis/R Code/Studie_FixRan/FixRanStudyGit.git/Imagen/MNI152.nii', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

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
### Calculate overlap
###############
##


# Overlap in each PAIR
OverlapP <- list(
	array(0,dim=c(numpoolmeths,NPAIRS),dimnames=list(poolmeth,paste('P',c(1:NPAIRS),sep=''))),
	array(0,dim=c(numpoolmeths,NPAIRS),dimnames=list(poolmeth,paste('P',c(1:NPAIRS),sep=''))),
	array(0,dim=c(numpoolmeths,NPAIRS),dimnames=list(poolmeth,paste('P',c(1:NPAIRS),sep=''))),
	array(0,dim=c(numpoolmeths,NPAIRS),dimnames=list(poolmeth,paste('P',c(1:NPAIRS),sep='')))
	)
names(OverlapP) <- metaMethods

# For loop over the PAIRS
for(p in 1:NPAIRS){
	# At run:
	print(paste('At pair ',p,sep=''))
	# Now we loop over the pooling and meta-analysis methods
	for(j in 1:numCols){
		if(grepl('ALE', METHODS[['MetaAnalysis']][j])){
			# Assign the MASKALE to the variable OVERLAPMASK: used when calculating overlap
			assign('OVERLAPMASK', MASKALE)
			# Load in the data: ALE uncorrected using matlab code of Eickhoff, then cluster corrected
			if(METHODS[['MetaAnalysis']][j] == "ALEUn"){
				# Group A: Z-values then going to P-values
				GroupA.Z <- readNIfTI(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/ALE/',METHODS[[2]][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				GroupA.tmp <- 1-pnorm(GroupA.Z, mean = 0, sd = 1)

				# Group B
				GroupB.Z <- readNIfTI(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/ALE/',METHODS[[2]][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				GroupB.tmp <- 1-pnorm(GroupB.Z, mean = 0, sd = 1)

				# Thresholding. Maps are already thresholded, hence NA to argument.
				GroupA <- ThreshPVal(GroupA.tmp, 0.001, mask = MASKALE)
				GroupB <- ThreshPVal(GroupB.tmp, 0.001, mask = MASKALE)
					rm(GroupA.tmp,GroupB.tmp,GroupA.Z,GroupB.Z)
			}else{
				# Group A: cFWE
				filePathA <- paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/ALE/',METHODS[['Pooling']][j],'/ALE/Results/',sep='')
				GroupA.tmp <- readNIfTI(dir(path = filePathA, pattern = '^OLS_cFWE05_001_', full = TRUE), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				GroupA <- ThreshALECluster(Thr_ALEMAP = GroupA.tmp)
					rm(GroupA.tmp, filePathA)

				# Group B: cFWE
				filePathB <- paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/ALE/',METHODS[['Pooling']][j],'/ALE/Results/',sep='')
				GroupB.tmp <- readNIfTI(dir(path = filePathB, pattern = '^OLS_cFWE05_001_', full = TRUE), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				GroupB <- ThreshALECluster(Thr_ALEMAP = GroupB.tmp)
					rm(GroupB.tmp, filePathB)
			}
		}else{
			# Assign the MASK to the variable OVERLAPMASK
			assign('OVERLAPMASK', MASK)
			# Load in PVal of group A: fixed and random effects MA
			load(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/PVal',sep=''))

			# Second approach:
			# Significance testing at P < 0.005 as well as Z > 1
				# Load in the Z-values of the MA
				# We need the placeholder name for either fixed or random effects MA, this is the first 3 letters of METHODS[[1]][j] (Fix or Ran)
				PLACEHOLDER <- METHODS[[1]][j] %>% substr(1,3)
				load(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/Zstat',PLACEHOLDER,sep=''))
			# The Z-values
			GroupA.Z <- get(paste('Zstat',PLACEHOLDER,sep=''))
				# Make identifier for Z > 1
				GroupAID.Z <- GroupA.Z > 1
			# Threshold P-values at 0.005, while using correct mask
			GroupA.P <- ThreshPVal(PVal, threshold = 0.005, mask = MASK)
			# Only select those with GroupAID.Z = TRUE
			GroupA.tmp <- GroupA.P
				GroupA.tmp[!GroupAID.Z] <- 0
			GroupA <- GroupA.tmp
				rm(PVal, GroupA.tmp, GroupAID.Z, PLACEHOLDER, GroupA.P)


			# Load in PVal of group B
			load(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/PVal',sep=''))

			# Second approach:
			# Significance testing at P < 0.005 as well as Z > 1
				# Load in the Z-values of the MA
				# We need the placeholder name for either fixed or random effects MA, this is the first 3 letters of METHODS[[1]][j] (Fix or Ran)
				PLACEHOLDER <- METHODS[[1]][j] %>% substr(1,3)
				load(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/Zstat',PLACEHOLDER,sep=''))
			# The Z-values
			GroupB.Z <- get(paste('Zstat',PLACEHOLDER,sep=''))
				# Make identifier for Z > 1
				GroupBID.Z <- GroupB.Z > 1
			# Threshold P-values at 0.005, while using correct mask
			GroupB.P <- ThreshPVal(PVal, threshold = 0.005, mask = MASK)
			# Only select those with GroupBID.Z = TRUE
			GroupB.tmp <- GroupB.P
				GroupB.tmp[!GroupBID.Z] <- 0
			GroupB <- GroupB.tmp
				rm(PVal, GroupB.tmp, GroupBID.Z, PLACEHOLDER, GroupB.P)

			}
		# Now that we have the thresholded map for group A and B, let us calculate the pairwise overlap
		method <- METHODS[[1]][j]
		pooling <- METHODS[[2]][j]
		OverlapP[[method]][pooling,p] <- OVERLAP(GroupA, GroupB, mask = OVERLAPMASK)
		}
	# Remove masks
	rm(OVERLAPMASK)
}

# Output the average overlap values over all iterations AND group level models
AvgOverlapOverRunsModels <- lapply(OverlapP, FUN = mean)
lapply(AvgOverlapOverRunsModels, FUN = round, digits = 2)

# Output the average overlap values over all iterations
AvgOverlapOverRuns <- lapply(OverlapP, FUN = rowMeans)
lapply(AvgOverlapOverRuns, FUN = round, digits = 2)

##
###############
### Plotting
###############
##

# Generic matrix
GM <- matrix(0, nrow=NRUNS,ncol=NRUNS)
	diag(GM) <- 1
	colnames(GM) <- rownames(GM) <- paste('l', 1:NRUNS, sep='')

# Correlations
CORRM <- c()

# Wrangle the correlation matrices
for(j in 1:numMetaMethods){
	for(i in 1:numpoolmeths){
		# Put the overlap values in a data frame
		index <- ((j - 1) * numpoolmeths) + i
		dat.tmp <- OverlapP[[j]][i,]

		# Transform to correlation matrix
		tmp.M <- GM
		tmp.M[lower.tri(tmp.M)] <- tmp.M[upper.tri(tmp.M)] <- dat.tmp
		tmp.M[lower.tri(tmp.M)] <- NA

		# Melt the data
		melted_M <- melt(tmp.M, na.rm = TRUE)

		# Add pooling method to it
		melted_M$pooling <- METHODS[[2]][index]

		# Add meta-analysis to it
		melted_M$MA <- METHODS[[1]][index]

		# Reset rownames
		rownames(melted_M) <- NULL

		# Add them together
		CORRM <- rbind(CORRM, melted_M)
	}
}

# Add labeling for plotting to factors of CORRM
CORRM$pooling <- factor(CORRM$pooling, levels = sort(poolmeth), labels = sort(ArtPOOLlabels))
CORRM$MA <- factor(CORRM$MA, levels = sort(metaMethods), labels = sort(ArtLABMA))


# Colour set: choose one!
col1 <- c('#08589e', '#7bccc4', '#ccebc5' )
col1 <- c('#7fc97f','#beaed4','#fdc086')
OverlapFigNumber <- c(9,"10A","10B")

# Overlap heatmap
quartz(height = 6.5, width = 6.5,
	type = 'tif', file = paste(LocFileSave, '/figure',OverlapFigNumber[WD],'_overlap_K',NSTUD,'.tif', sep = ''), bg = 'white', canvas = 'white')
ggplot(data = CORRM, aes(Var2, Var1, fill = value))+
geom_tile(color = "white") +
facet_grid(pooling ~ MA) +
geom_text(aes(label = str_replace(as.character(round(value,2)), "^0\\.", ".")), colour = "white", size = ifelse(WD == 1,2.5,4)) +
scale_fill_gradient2(low = col1[1], mid = col1[2], high = col1[3],
  midpoint = 0.5, limit = c(0,1), space = "Lab",
  name="Overlap\nCoefficient") +
	scale_x_discrete(name=paste("ITERATION (K = ", NSTUD, ")", sep=""),
			labels = 1:NRUNS) +
     scale_y_discrete(name="ITERATION",
		 	labels = 1:NRUNS) +
 theme_minimal() +
theme(axis.text.x = element_text(angle = 0, vjust = 1,
   		size = 10, hjust = 1),
   legend.box.margin = margin(0, 10, 0, 10),
	 legend.position = 'top',
	 legend.justification = "left",
	 legend.title = element_text(size = 8.5),
	 legend.text.align = 0.40,
	 plot.margin = unit(c(0,0,0,0), 'mm'),
	 panel.border = element_blank())+
coord_fixed()
dev.off()




##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
