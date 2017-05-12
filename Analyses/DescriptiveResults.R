####################
#### TITLE:     Calculate descriptive results for activation reliability of the MAs.
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
WD <- 1

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


# Mask for thresholding the PMaps of fixed and random effects MA.
MASK <- readNIfTI(paste(WDs[[WD]],'/Run_1/GroupAnalysis/mask.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
# Mask for ALE
MASKALE <- readNIfTI('MNI152.nii', verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

# Labels to be used in paper for group level models
ArtPOOLlabels <- c('OLS', 'Fixed Effects', 'Mixed Effects')

# Labels to be used in paper for meta-analyses
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
### Descriptive results
###############
##

# Apply 26 point searching cluster algorithm of FSL to the thresholded maps.
# Looking at amount of clusters, size of them,  amount of overlapping clusters, ...

# Extra library
library(tables)

# Function to have custom cluster command
CustomCluster <- function(BrainMap, TempWD, OutputName){
	# First save the current wd (to get there back later on)
	currentWD <- getwd()

	# Name of the Brain Map
	NameMap <- deparse(substitute(BrainMap))

	# Write the BrainMap to a nifti file in a temporary directory
	writeNIfTI(BrainMap, file = paste(TempWD,'/',NameMap, sep=''), gzipped = FALSE)

	# Create command for clustering in FSL
	command <- paste('/usr/local/fsl/bin/cluster -i ', NameMap,'.nii -t 1 -o ', OutputName, sep='')
	# Change wd to the tempWD folder, put system environment to NIFTI and then execute command
	setwd(TempWD)
	Sys.setenv(FSLOUTPUTTYPE="NIFTI")
	system(command)

	# Read in this image to return it
	ResultImage <- readNIfTI(paste(TempWD,'/', OutputName, '.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]

	# Create command to clean temp folder
	commandRM <- paste('rm ', TempWD, '/',OutputName,'.nii ', TempWD,'/',NameMap,'.nii', sep = '')
	system(commandRM)

	# Go back to original WD
	setwd(currentWD)

	# Return
	return(ResultImage)
}

# Object with information
ClustInfo <- c()


# For loop over the PAIRS
for(p in 1:NPAIRS){
	# At run:
	print('------------------------------------------------')
	print(paste('AT PAIR ',p,sep=''))
	print('------------------------------------------------')
	# Now we loop over the pooling and meta-analysis methods
	for(j in 1:numCols){
		if(grepl('ALE', METHODS[['MetaAnalysis']][j])){
			# Assign ALE dimension to the variable SwitchDIM
			SwitchDIM <- dim(MASKALE)
			# Load in the data: ALE uncorrected using matlab code of Eickhoff, then cluster corrected
			if(METHODS[['MetaAnalysis']][j] == "ALEUn"){
				# Group A: read in Z-values and transform to P-values
				GroupA.Z <- readNIfTI(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/ALE/',METHODS[[2]][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				GroupA.tmp <- 1-pnorm(GroupA.Z, mean = 0, sd = 1)

				# Group B
				GroupB.Z <- readNIfTI(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/ALE/',METHODS[[2]][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				GroupB.tmp <- 1-pnorm(GroupB.Z, mean = 0, sd = 1)

				# Thresholding and masking
				# To have descriptive results, we use the mask specific to each run!
					# As opposed to intersection in overlap calculations.
					MASKA <- MASKB <- array(1, dim=SwitchDIM)
					MASKA[is.na(GroupA.Z)] <- 0; MASKB[is.na(GroupB.Z)] <- 0
				GroupA <- ThreshPVal(GroupA.tmp, 0.001, mask = MASKA)
				GroupB <- ThreshPVal(GroupB.tmp, 0.001, mask = MASKB)
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
				# Assign IMAGEN dimension to the variable SwitchDIM
				SwitchDIM <- DIM


				# Load in PVal of group A: fixed and random effects MA
				load(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/PVal',sep=''))

				# Significance testing at P < 0.005 as well as Z > 1
					# Load in the Z-values of the MA
					# We need the placeholder name for either fixed or random effects MA, this is the first 3 letters of METHODS[[1]][j] (Fix or Ran)
					PLACEHOLDER <- METHODS[[1]][j] %>% substr(1,3)
					load(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/Zstat',PLACEHOLDER,sep=''))
				# The Z-values
				GroupA.Z <- get(paste('Zstat',PLACEHOLDER,sep=''))
					# Make identifier for Z > 1
					GroupAID.Z <- GroupA.Z > 1
				# Threshold P-values at 0.005, while using mask: universal mask of IMAGEN
				GroupA.P <- ThreshPVal(PVal, threshold = 0.005, mask = MASK)
				# Only select those with GroupAID.Z = TRUE
				GroupA.tmp <- GroupA.P
					GroupA.tmp[!GroupAID.Z] <- 0
				GroupA <- GroupA.tmp
					rm(PVal, GroupA.tmp, GroupAID.Z, PLACEHOLDER, GroupA.P)


				# Load in PVal of group B
				load(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/PVal',sep=''))

				# Significance testing at P < 0.005 as well as Z > 1
					# Load in the Z-values of the MA
					# We need the placeholder name for either fixed or random effects MA, this is the first 3 letters of METHODS[[1]][j] (Fix or Ran)
					PLACEHOLDER <- METHODS[[1]][j] %>% substr(1,3)
					load(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/Zstat',PLACEHOLDER,sep=''))
				# The Z-values
				GroupB.Z <- get(paste('Zstat',PLACEHOLDER,sep=''))
					# Make identifier for Z > 1
					GroupBID.Z <- GroupB.Z > 1
				# Threshold P-values at 0.005, while using mask: universal mask of IMAGEN
				GroupB.P <- ThreshPVal(PVal, threshold = 0.005, mask = MASK)
				# Only select those with GroupBID.Z = TRUE
				GroupB.tmp <- GroupB.P
					GroupB.tmp[!GroupBID.Z] <- 0
				GroupB <- GroupB.tmp
					rm(PVal, GroupB.tmp, GroupBID.Z, PLACEHOLDER, GroupB.P)

			}

		# Now that we have the thresholded maps, let us see how clusters of voxels overlap
		method <- METHODS[[1]][j]
		pooling <- METHODS[[2]][j]

		# Clustering GroupA
		GroupAClust <- CustomCluster(BrainMap = GroupA, TempWD = paste(WDs[[WD]], '/temp', sep=''), OutputName = 'GroupAClust')

		# Clustering GroupB
		GroupBClust <- CustomCluster(BrainMap = GroupB, TempWD = paste(WDs[[WD]], '/temp', sep=''), OutputName = 'GroupBClust')

		# Amount of clusters in each image
		NumClustA <- max(GroupAClust)
		NumClustB <- max(GroupBClust)

		# Average voxel size of the clusters
		AvVoxA <- sum(GroupAClust != 0)/NumClustA
		AvVoxB <- sum(GroupBClust != 0)/NumClustB

		# Overlapping clusters
		ArrayA <- array(GroupAClust, dim = prod(SwitchDIM))
		ArrayB <- array(GroupBClust, dim = prod(SwitchDIM))
			# Using for loop to check over every voxel whether there is a value in
			# one group or the other. If so, then we delete this cluster.
			# Remaining clusters are hence unique clusters (that do not overlap).
			for(i in 1:prod(SwitchDIM)){
				NumA <- ArrayA[i]
				NumB <- ArrayB[i]
				SumAB <- NumA + NumB
				if(SumAB > NumA) ArrayA[ArrayA==NumA] <- 0
				if(SumAB > NumB) ArrayB[ArrayB==NumB] <- 0
			}

		# Determine the overlapping clusters
		IDAoverlap <- array(GroupAClust, dim = prod(SwitchDIM)) != ArrayA
		IDBoverlap <- array(GroupBClust, dim = prod(SwitchDIM)) != ArrayB

		# Number of overlapping clusters in each group (excluding zeroes)
		OverlapClustA <- length(unique(array(GroupAClust, dim = prod(SwitchDIM))[IDAoverlap]))
		OverlapClustB <- length(unique(array(GroupBClust, dim = prod(SwitchDIM))[IDBoverlap]))

		# Average cluster size of overlapping clusters in each group
		OverlapClustASize <- length(array(GroupAClust, dim = prod(SwitchDIM))[IDAoverlap]) / OverlapClustA
		OverlapClustBSize <- length(array(GroupBClust, dim = prod(SwitchDIM))[IDBoverlap]) / OverlapClustB

		# Number of unique clusters in each group (excluding the zeroes)
		UniqueClustA <- length(unique(ArrayA[ArrayA != 0]))
		UniqueClustB <- length(unique(ArrayB[ArrayB != 0]))

		# Cluster size of the unique clusters
		UniqueClustASize <- sum(ArrayA != 0)/UniqueClustA
		UniqueClustBSize <- sum(ArrayB != 0)/UniqueClustB

		# Amount of non overlapping clusters with > average cluster size
		LargClustA <- sum(table(ArrayA[ArrayA != 0]) > AvVoxA)
		LargClustB <- sum(table(ArrayB[ArrayB != 0]) > AvVoxB)


		# Save all the numbers in a temp data frame
		ClustInfo.tmp <- data.frame('Amount Of Clusters' = c(NumClustA, NumClustB),
																'Average Voxel Size in Cluster' = c(AvVoxA, AvVoxB),
																'Overlapping Clusters' = c(OverlapClustA,OverlapClustB),
																'Avg Vox Size Overlap' = c(OverlapClustASize,OverlapClustBSize),
																'Unique clusters' = c(UniqueClustA, UniqueClustB),
																'Avg Vox Size Unique' = c(UniqueClustASize, UniqueClustBSize),
																'Large unique clusters' = c(LargClustA, LargClustB),
																'Group' = c('A','B'),
																'Method' = method,
																'Pooling' = pooling,
																'Run' = c(PAIRS[p,1],PAIRS[p,2]))
		# Combine in one data frame
		ClustInfo <- rbind(ClustInfo, ClustInfo.tmp)
		}
	# Remove masks
	rm(MASKA, MASKB, SwitchDIM)
}

# Save or load object
#
#save(ClustInfo, file = paste(WDs[[WD]], '/ClustInfo.RData', sep =''))
load(file=paste(WDs[[WD]], '/ClustInfo.RData', sep =''))


# View Cluster Info
ClustInfo
str(ClustInfo)
# Let us create useful aggregated tables: CARE AT THIS STEP!
		# Consider the case in which we caluclated 7 iterations of the design.
		# We only have 7 resulting images (per pooling method and MA) for which we can caluclate the amount of clusters and voxels per cluster.
		# However(!).
		# The table we get follows the PAIRS object (which is pairwise combination of the iterations)
		# Hence we have unbalanced duplicates (e.g. we compare iteration 1 with all seven others so we have 6 times iteration 1 in the table, but iteration 6 only with 7 once).
		# These pairwise comparissons are useful for calculating amount of unique clusters and uniqe large clusters.
		# We thus first need to calculate the average of unique clusters and unique large clusters.
		# Then can we select the unique iterations and calculate the average amount of clusters and average voxel size!
# First aggregate while keeping group A and B (to check if they are more or less similar)
AggClustInfoAB <- aggregate(cbind(Unique.clusters, Large.unique.clusters) ~ Group + Method + Pooling, data = ClustInfo, FUN=mean)

# OK. Now over both groups as well
AggClustInfoAgg <- aggregate(cbind(Unique.clusters, Large.unique.clusters) ~ Method + Pooling, data = ClustInfo, FUN=mean)
AggClustInfo <- tabular((Pooling * Method) ~ (n=1) + Format(digits=2) * (Unique.clusters + Large.unique.clusters) * (mean + sd), data = ClustInfo)

# Now we can select the unique iterations (which correspond to NRUNS)
UniqueIDs <- rownames(ClustInfo) %in% rownames(unique(ClustInfo[,c('Method','Pooling','Run')]))
UniqueClustInfo <- ClustInfo[UniqueIDs, ]
# And aggregate
AggUniqeClustAgg <- aggregate(cbind(Amount.Of.Clusters, Average.Voxel.Size.in.Cluster) ~ Method + Pooling, data = UniqueClustInfo, FUN=mean)
AggUniqeClust <- tabular((Pooling * Method) ~ (n=1) + Format(digits=2) * (Amount.Of.Clusters + Average.Voxel.Size.in.Cluster) * (mean + sd), data = UniqueClustInfo)

# Combine both tables
ClustTable <- cbind(AggUniqeClust, AggClustInfo[,c(2:5)])
ClustTable
	# Using second approach
	ClustTableAgg <- cbind(AggUniqeClustAgg, AggClustInfoAgg[,c(3:4)])
	# Mean over MA
	aggregate(cbind(Amount.Of.Clusters, Average.Voxel.Size.in.Cluster, Unique.clusters, Large.unique.clusters) ~ Pooling, data = ClustTableAgg, FUN = mean)


##########################################
## Let us combine the cases into one table
##########################################

# Renaming for printing purpose
labelingPooling <- c('Fixed Effects', 'OLS', 'Mixed Effects')
labelingMA <- c('Fixed Effects MA', 'Random Effects MA','ALE cFWE','ALE Uncorrected')

# Working directory indices, number of studies and empyt vector
	wdsClustInfo <- c(6,5,7)				# this order as this is K=10, 20, 35
	KforClustInfo <- c(10, 20, 35)
	AllClustTable <- c()

# Using for loop: loading and then adding the data
for(i in 1:length(wdsClustInfo)){
	# Load in ClustInfo
	load(file=paste(WDs[[wdsClustInfo[i]]], '/ClustInfo.RData', sep =''))

	# Add K studies column
	ClustInfo$KSTUD <- KforClustInfo[i]

	# Renaming for printing purpose and make K a factor
	ClustInfo$Pooling <- factor(ClustInfo$Pooling, labels = labelingPooling)
	ClustInfo$Method <- factor(ClustInfo$Method, labels = labelingMA[c(1,2,4,3)])
	ClustInfo$KSTUD <- factor(ClustInfo$KSTUD)

	# Aggregate unique clusters and unique large clusters
	AggClustInfo <- tabular(( RowFactor(KSTUD, spacing = 1) * RowFactor(Pooling, spacing = 1, space = 0.5) * Method) ~ (n=1) + Format(digits=2) * (Unique.clusters + Large.unique.clusters) * (mean + sd), data = ClustInfo)

	# Select the uniqe iterations
	UniqueIDs <- rownames(ClustInfo) %in% rownames(unique(ClustInfo[,c('Method','Pooling','Run')]))
	UniqueClustInfo <- ClustInfo[UniqueIDs, ]

	# Aggregate the amount of clusters and the average voxel size per cluster
	AggUniqeClustInfo <- tabular(( RowFactor(KSTUD, spacing = 1) * RowFactor(Pooling, spacing = 1, space = 0.5) * Method) ~ (n=1) + Format(digits=2) * (Amount.Of.Clusters + Average.Voxel.Size.in.Cluster) * (mean + sd), data = UniqueClustInfo)

	# Combine AggClustInfo and AggUniqeClustInfo
	ClustTable <- cbind(AggUniqeClustInfo, AggClustInfo[,c(2:5)])

	# It is kind of tricky to add calculations to the tables package.
	# Hence I will just output in .txt file the percentages and then add them manually to latex
	percentages <- as.numeric(as.matrix(ClustTable[,8])[c(3:14),4]) /
					as.numeric(as.matrix(ClustTable[,2])[c(3:14),4])
		write.table(round(percentages, 2), file = paste(WDs[[wdsClustInfo[i]]], '/percentages_', KforClustInfo[i],'.txt', sep =''),
							col.names = FALSE, row.names = FALSE)

	# Add data frames of K = 10, 20 and 35
	AllClustTable <- rbind(AllClustTable, ClustTable)
	rm(ClustInfo, AggClustInfo, UniqueIDs, UniqueClustInfo, AggUniqeClustInfo, ClustTable)
}
latex(AllClustTable)


# Print (1) average over all pooling methods and (2) average over all MA's
TableWithAllClusterInfo <- c()
for(i in 1:length(wdsClustInfo)){
	# Load in ClustInfo
	load(file=paste(WDs[[wdsClustInfo[i]]], '/ClustInfo.RData', sep =''))
	# Add K studies column
	ClustInfo$KSTUD <- KforClustInfo[i]
	# Join the data
	TableWithAllClusterInfo <- rbind(TableWithAllClusterInfo, ClustInfo)
}
# Average over group level models
print('Average over group level models')
aggregate(cbind(Amount.Of.Clusters, Average.Voxel.Size.in.Cluster, Unique.clusters, Large.unique.clusters) ~ Method + KSTUD, data = TableWithAllClusterInfo, FUN = mean)

# Average over MA's
print('Average over MAs')
aggregate(cbind(Amount.Of.Clusters, Average.Voxel.Size.in.Cluster, Unique.clusters, Large.unique.clusters) ~ Pooling + KSTUD, data = TableWithAllClusterInfo, FUN = mean)




##########################################
## Make plots of these numbers: APPENDIX
##########################################

ClustInfoTbl <- tbl_df(ClustInfo)
# Just a check for the amount of clusters
# ClustInfoTbl %>% mutate(test = Overlapping.Clusters + Unique.clusters,
# 												check = Amount.Of.Clusters == test) %>% select(check) %>% summarise(sum(check))


# Test case: pick your method
SelectedInfo <- ClustInfoTbl %>% filter(Method == 'FixedEffUn' & Pooling == 'FixedEffect') %>%
									select(Overlapping.Clusters, Unique.clusters,Avg.Vox.Size.Overlap, Avg.Vox.Size.Unique)

LongSelectedInfo <- data.frame('AmountClust' = c(SelectedInfo$Overlapping.Clusters, SelectedInfo$Unique.clusters),
					'ClustSize' = c(SelectedInfo$Avg.Vox.Size.Overlap, SelectedInfo$Avg.Vox.Size.Unique),
					'Source' = rep(c('Overlapping', 'Unique'), each = length(SelectedInfo$Overlapping.Clusters))) %>%
					tbl_df()

ggplot(LongSelectedInfo, aes(x = AmountClust)) + 	geom_histogram(data = LongSelectedInfo %>% select(AmountClust), fill = 'grey', bins = 25) +
	geom_histogram(bins = 25) + facet_grid(Source ~ .) +
		scale_x_continuous(breaks = waiver(), labels = waiver(), name = 'Amount of clusters')


# Now do them all!
# Make data frame with all pooling and MA methods
SelectedInfo <- ClustInfoTbl %>% select(Overlapping.Clusters, Unique.clusters,Avg.Vox.Size.Overlap, Avg.Vox.Size.Unique, Method, Pooling)

# We will need to stack the overlapping values on the unique values in one vector (for plotting)
# Hence we copy two times the method and pooling name
MethodSelect <- bind_rows('Method' = data.frame(SelectedInfo$Method), 'Method' = data.frame(SelectedInfo$Method)) %>%
									rename(Method = SelectedInfo.Method) %>% tbl_df()
MethodSelect$Method <- factor(MethodSelect$Method, levels = metaMethods, labels = ArtLABMA[c(1,2,4,3)])

PoolingSelect <- bind_rows('Pooling' = data.frame(SelectedInfo$Pooling), 'Pooling' = data.frame(SelectedInfo$Pooling)) %>%
									rename(Pooling = SelectedInfo.Pooling) %>% tbl_df()
PoolingSelect$Pooling <- factor(PoolingSelect$Pooling, levels = poolmeth, labels = ArtPOOLlabels[c(2,1,3)])

# Now we make a data frame in long format (overlapping and unique values) onto each other.
LongSelectedInfo <- data.frame('AmountClust' = c(SelectedInfo$Overlapping.Clusters, SelectedInfo$Unique.clusters),
					'ClustSize' = c(SelectedInfo$Avg.Vox.Size.Overlap, SelectedInfo$Avg.Vox.Size.Unique),
					'Source' = rep(c('Overlapping', 'Unique'), each = length(SelectedInfo$Overlapping.Clusters)),
					'Method' = MethodSelect,
					'Pooling' = PoolingSelect) %>%
					mutate(CombName = paste(as.character(Pooling), as.character(Method), sep = ' - ')) %>%
					tbl_df()

# Histogram over all data
ggplot(data = LongSelectedInfo %>% select(AmountClust, Source),
			aes(x = AmountClust)) + geom_histogram(bins = 25, fill = 'grey') + facet_grid(Source ~ .)

# Now we can make for each of the combination group and MA a plot (using do(plots ))
PlotsSelAmount <- LongSelectedInfo %>% group_by(Method, Pooling) %>% do(plots = ggplot(data = .) +
										aes(x = AmountClust) +
										geom_histogram(data = LongSelectedInfo %>% select(AmountClust, Source), fill = 'grey', bins = 25) +
											geom_histogram(bins = 25) + facet_grid(Source ~ .) +
												scale_x_continuous(breaks = waiver(), labels = waiver(), name = 'Amount of [Uniq./Over.] clusters per map') +
												ggtitle(unique(.$CombName)))
# To plot them, use
PlotsSelAmount$plots

# Same for cluster sizes!
PlotsSelSize <- LongSelectedInfo %>% group_by(Method, Pooling) %>% do(plots = ggplot(data = .) +
										aes(x = ClustSize) +
										geom_histogram(data = LongSelectedInfo %>% select(ClustSize, Source), fill = 'grey', bins = 25) +
										geom_histogram(bins = 25) + facet_grid(Source ~ .) +
												scale_x_continuous(breaks = waiver(), labels = waiver(), name = 'Size of [Uniq./Over.] clusters (in voxels)') +
												ggtitle(unique(.$CombName)))
# To plot them, use
PlotsSelSize$plots

# Saving
# 2 x 2
for(i in 1:length(PlotsSelSize$plots)){
	quartz(height = 6, width = 10, type = 'png',
		file = paste(LocFileSave, '/DescriptivePlots/', NSTUD, '/Des_',i,'K_',NSTUD, '.png', sep = ''),
		dpi = 600, bg = 'white', canvas = 'white')
	grid.arrange(PlotsSelAmount$plots[[i]], PlotsSelSize$plots[[i]], ncol = 2)
	dev.off()
}

##################
# Without facetting according to overlapping/unique
##################

PlotsSelAmountOverlaid <- LongSelectedInfo %>% group_by(Pooling) %>% rename(Type = Source) %>% do(plots = ggplot(data = .) +
										aes(x = AmountClust) +
											geom_histogram(bins = 25, alpha = 0.5, position = 'identity', aes(fill = Type)) + facet_grid(Method ~ .) +
												scale_x_continuous(breaks = waiver(), labels = waiver(), name = 'Amount of [Uniq./Over.] clusters per map') +
												ggtitle(paste0('Second level model: ', unique(.$Pooling))))

PlotsSelAmountOverlaid$plots


PlotsSelSizeOverlaid <- LongSelectedInfo %>% group_by(Pooling) %>% rename(Type = Source) %>% do(plots = ggplot(data = .) +
										aes(x = ClustSize) +
											geom_histogram(bins = 25, alpha = 0.5, position = 'identity', aes(fill = Type)) + facet_grid(Method ~ .) +
												scale_x_continuous(breaks = waiver(), labels = waiver(), name = 'Size of [Uniq./Over.] clusters (in voxels)') +
												ggtitle(paste0('Second level model: ', unique(.$Pooling))))

PlotsSelSizeOverlaid$plots

# Saving
# 2 x 2
for(i in 1:length(PlotsSelSizeOverlaid$plots)){
	quartz(height = 6, width = 10, type = 'png',
		file = paste(LocFileSave, '/DescriptivePlots/', NSTUD, '/OverLaidDes_',i,'_K_',NSTUD, '.png', sep = ''),
		dpi = 600, bg = 'white', canvas = 'white')
	grid.arrange(PlotsSelAmountOverlaid$plots[[i]], PlotsSelSizeOverlaid$plots[[i]], ncol = 2)
	dev.off()
}




