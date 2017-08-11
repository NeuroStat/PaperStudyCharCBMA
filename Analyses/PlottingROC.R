####################
#### TITLE:     Calculate ROC curves when comparing MA's with the reference.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/Analyses
#### First Modified: 03/05/2016
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
library(ggplot2)
library(ggthemes)
library(oro.nifti)
library(dplyr)
library(tidyr)

# Function to caluclate ROC values based on P-map
ROC_PMAP <- function(GTmap, PValmap,number.thresholds=100){
	# The thresholds, the GT in one array and the p-map in one array
	thresholds <- seq(0,1,length.out=number.thresholds)
	#thresholds <- seq(0,1,by=0.01)
	GT <- array(GTmap, dim = prod(dim(GTmap)))
	pmap <- array(PValmap, dim = prod(dim(PValmap)))

	# For loop over the thresholds to calculate false and true positives
	TFP <- TTP <- c()
	for(t in 1:length(thresholds)){
		# First (re)-threshold the maps
		Tpmap <- array(0, dim = prod(dim(PValmap)))
		idTpmap <- pmap <= thresholds[t]
		Tpmap[idTpmap] <- 1

		# False positive rate: if GT = 0, and p-map = 1, then difference equals -1
		FP <- round(sum((GT - Tpmap) == -1, na.rm = TRUE) / sum(GT == 0, na.rm=TRUE), 3)
		# True positive rate: if both GT and p-map have 1, then sum = 2
		TP <- round(sum((GT + Tpmap) == 2, na.rm = TRUE) / sum(GT, na.rm = TRUE), 3)

		# Save in vector (thresholdFP, thresholdTP)
		TFP <- c(TFP, FP)
		TTP <- c(TTP, TP)
	}
	# Save in data.frame
	ROC <- data.frame('value' = c(TFP,TTP), 'type' = rep(c('FP', 'TP'), each = length(thresholds)), 'threshold' = round(rep(thresholds, 2),2))
	# Return it
	return(ROC)
}

# Function to caluclate area under curve
comp_auc <- function(FP,TP){
	# Calculate area under curve.
		# x-axis = FPR
		# y-axis = TPR
	auc <- 0
	for (i in 2:length(FP)) {
		auc <- auc + 0.5 * (FP[i] - FP[i-1]) * (TP[i] + TP[i-1])
	}
	return(data.frame(AUC = auc))
}

# Function to calculate either the pID (assumption of independence or positive dependence (hence no negative correlation))
#     and the pD which makes no assumption about the correlation of the data
FDR <- function(p,q,mask){
  if(!all(dim(p) == dim(mask))){stop("Dimension of p-value image is not equal to mask!")}
  if(!all(mask %in% c(0,1))){stop("Mask can only contain 1 and 0 values!")}
  p <- array(p,dim=prod(dim(p)))
  mask <- array(mask,dim=prod(dim(mask)))
  p[mask==0] <- NA
  p <- p[!is.na(p)]
  p <- sort(p)
  V <- length(p)
  I <- c(1:V)

  cVID <- 1
  cVN <- sum(1/(1:V))

  pID <- p<=((I/V)*(q/cVID))
  pN <- p<=((I/V)*(q/cVN))

  ThpID <- max(p[pID])
  ThpN <- max(p[pN])

  thresholds <- list(ThpID,ThpN)
  names(thresholds) <- c('pID','pN')
  return(thresholds)

}

# Number of runs.
if(WD == 1) NRUNS <- 7
if(WD == 2) NRUNS <- 3
if(WD == 3) NRUNS <- 2


# Number of studies in each meta-analysis
NSTUD_FU <- function(WD){
	if(WD == 1) NSTUD <- 10
	if(WD == 2) NSTUD <- 20
	if(WD == 3) NSTUD <- 35
	return(NSTUD)
}
NSTUD <- NSTUD_FU(WD)


# Dimension of the data
DIM <- c(53,63,46)

# Vector of pooling methods
poolmeth <- c(
	'FixedEffect',
	'OLS',
	'MixedEffect'
	)
numpoolmeths <- length(poolmeth)

# Labels to be used in paper for group level models
ArtPOOLlabels <- c('OLS', 'Fixed Effects', 'Mixed Effects')

# Labels of meta-analyses
ArtLABMA <- c('ALE', 'Fixed Effects MA','Random Effects MA')


##
###############
### Calculate numbers
###############
##

# Array of meta-analysis methods
metaMethods <- c('FixedEffUn', 'RanEffUn','ALE')
numMetaMethods <- length(metaMethods)

# number of columns in the data set
numCols <- numMetaMethods * numpoolmeths

# number of thresholds considered in the ROC
number.thresholds <- 100

# List of all the pooling methods and meta-analyses when loading data
METHODS <- list(
	array(rep(metaMethods,each=numpoolmeths)),
	array(rep(poolmeth,numMetaMethods)))
names(METHODS) <- c('MetaAnalysis', 'Pooling')

# False Positive Rate and True Positive Rate
FPR <- TPR <- FPR_Un001 <- TPR_Un001 <- FPR_FDR05 <- TPR_FDR05 <- rep(
	list(
		array(0,dim=c(number.thresholds,NRUNS),dimnames=list(c(1:number.thresholds),paste('R',c(1:NRUNS),sep='')))
		)	,numCols
	)
	# Naming structure: CARE!!
			# The order = each MA first, then the pooling methods. So for fixed effect MA we see the 3 pooling methods, etc...
			# THIS NEED TO CORRESPOND WITH THE METHODS OBJECT!
	names_two <- unique(expand.grid(rev(METHODS)))
	names_one <- paste(names_two[,1],names_two[,2], sep=":")
names(FPR) <- names(TPR) <- names(FPR_Un001) <- names(TPR_Un001) <- names(FPR_FDR05) <- names(TPR_FDR05) <- names_one

# We will need to for loop over runs
for(r in  1:NRUNS){
	# At run:
	print(paste('At run ',r,sep=''))
	# Working directory for this run
	runWD <- paste(WDs[[WD]],'/Run_',r,'/MetaAnalyses',sep='')

	# Load and mask the GT: two versions, one for ALE and one for fixed and random effects MA
	###########
	GroupMapFR  <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/thresh_zstat1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
	GroupMaskFR <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/mask.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
	idGroupMaskFR <- GroupMaskFR==0
		GroupMapFR[idGroupMaskFR] <- NA
	###########
	GroupMapALE <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/flirt_IMAGEN_to_MNI_then_thresh_zstat1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
		GroupMapALE[IDMNI152] <- NA

	# Second version of GT: threshold at uncorrected < 0.001
	# Load and mask the GT: two versions, one for ALE and one for fixed and random effects MA
	###########
	UnZGroupMapFR  <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/stats/zstat1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
	GroupMaskFR <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/mask.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
	idGroupMaskFR <- GroupMaskFR==0
		UnZGroupMapFR[idGroupMaskFR] <- NA
	Un001GroupMapFR <- array(0, dim = DIM)
		Un001GroupMapFR[idGroupMaskFR] <- NA
	Un001GroupMapFR[UnZGroupMapFR >= qnorm(0.001, lower.tail = FALSE)] <- 1
	###########
	UnZGroupMapALE <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/stats/flirted_IMAGEN_to_MNI_zstat1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
		UnZGroupMapALE[IDMNI152] <- NA
	Un001GroupMapALE <- array(0, dim = dim(UnZGroupMapALE))
		Un001GroupMapALE[IDMNI152] <- NA
	Un001GroupMapALE[UnZGroupMapALE >= qnorm(0.001, lower.tail = FALSE)] <- 1

	# Third version of GT: threshold at FDR < 0.05
	# Load and mask the GT: two versions, one for ALE and one for fixed and random effects MA
	###########
	FDR05GroupMapFR <- array(0, dim = DIM)
		FDR05GroupMapFR[idGroupMaskFR] <- NA
	UNPGroupMapFR <- pnorm(UnZGroupMapFR, lower.tail = FALSE)
		# Calculate the FDR threshold, assuming positive correlation
		PThreshFR <- FDR(p = UNPGroupMapFR, q = 0.05, mask = GroupMaskFR)$pID
	FDR05GroupMapFR[UNPGroupMapFR <= PThreshFR] <- 1

	###########
	# UnZGroupMapALE <- readNIfTI(paste(WDs[[WD]],'/Run_',r,'/GroupAnalysis/stats/flirted_IMAGEN_to_MNI_zstat1.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
	# 	UnZGroupMapALE[IDMNI152] <- NA
	FDR05GroupMapALE <- array(0, dim = dim(UnZGroupMapALE))
		FDR05GroupMapALE[IDMNI152] <- NA
		# Get p-values
		UNPGroupMapALE <- pnorm(UnZGroupMapALE, lower.tail = FALSE)
			# Calculate the FDR threshold, assuming positive correlation
			BinaryALE_MASK <- array(0, dim = dim(UnZGroupMapALE))
				BinaryALE_MASK[!IDMNI152] <- 1
	PThreshALE <- FDR(p = UNPGroupMapALE, q = 0.05, mask = BinaryALE_MASK)$pID
	FDR05GroupMapALE[UNPGroupMapALE <= PThreshALE] <- 1

	# Load the P-map (and/or Z-map) of the considered meta-analyses in each pooling scenario.
	for(j in 1:numCols){
		# P-value: ALE has different naming structure!
		if(grepl('ALE', METHODS[['MetaAnalysis']][j])){
			# Assign GroupMapALE to GroupMap
			assign('GroupMap', GroupMapALE)
			assign('GroupMap_Un001', Un001GroupMapALE)
			assign('GroupMap_FDR05', FDR05GroupMapALE)
			# Load in the data: ALE uncorrected z-values: still need to transfrom these to P-values!
				ZVal.tmp <- readNIfTI(paste(runWD,'/',METHODS[[1]][j],'/',METHODS[[2]][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
				PVal.tmp <- 1-pnorm(ZVal.tmp, mean = 0, sd = 1)
					# Put NA values in PVal to NA in the GroupMap and to 1 in the PVal
					PVal <-	PVal.tmp
					GroupMap[is.na(PVal)] <- NA
					PVal[is.na(PVal)] <- 1
		}else{
			# Assign GroupMapFR to GroupMap
			assign('GroupMap', GroupMapFR)
			assign('GroupMap_Un001', Un001GroupMapFR)
			assign('GroupMap_FDR05', FDR05GroupMapFR)
			load(paste(runWD,'/',METHODS[[1]][j],'/',METHODS[[2]][j],'/PVal',sep='')); if(!exists('PVal')) print('WARNING, P-VALUE NOT FOUND')

		}

		# ROC values
		ROC_values <- ROC_PMAP(GroupMap, PVal,number.thresholds = number.thresholds)
		ROC_values_Un001 <- ROC_PMAP(GroupMap_Un001, PVal,number.thresholds = number.thresholds)
		ROC_values_FDR05 <- ROC_PMAP(GroupMap_FDR05, PVal,number.thresholds = number.thresholds)

		# Gather the FP's and TP's in separate matrices
		FPR[[names_one[j]]][,r] <- ROC_values[which(ROC_values$type == 'FP'),'value']
		TPR[[names_one[j]]][,r] <- ROC_values[which(ROC_values$type == 'TP'),'value']

			# Uncorrected 001
			FPR_Un001[[names_one[j]]][,r] <- ROC_values_Un001[which(ROC_values_Un001$type == 'FP'),'value']
			TPR_Un001[[names_one[j]]][,r] <- ROC_values_Un001[which(ROC_values_Un001$type == 'TP'),'value']

			# FDR 05
			FPR_FDR05[[names_one[j]]][,r] <- ROC_values_FDR05[which(ROC_values_FDR05$type == 'FP'),'value']
			TPR_FDR05[[names_one[j]]][,r] <- ROC_values_FDR05[which(ROC_values_FDR05$type == 'TP'),'value']

		# Remove pmap
		rm(PVal, ROC_values, GroupMap)
	}
}

# Take mean rowwise, but calculate SD as well
MEANFPR.tmp <- lapply(FPR,rowMeans)
	MEANFPR <- t(do.call(rbind, MEANFPR.tmp))
MEANTPR.tmp <- lapply(TPR,rowMeans)
SDTPR.tmp <- lapply(TPR, function(x) apply(x, 1, sd))
	MEANTPR <- t(do.call(rbind, MEANTPR.tmp))
	SDTPR <- t(do.call(rbind, SDTPR.tmp))

# Put them in data frames
FPR.DATA <- data.frame('value' = matrix(MEANFPR,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

TPR.DATA <- data.frame('value' = matrix(MEANTPR,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

SDTPR.DATA <- data.frame('value' = matrix(SDTPR,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

TFPR.DATA <- data.frame('FP' = matrix(MEANFPR,ncol=1),
						'TP' = matrix(MEANTPR,ncol=1),
						'SDTP' = matrix(SDTPR,ncol=1),
						'pooling' = array(rep(rep(poolmeth,each=number.thresholds),numMetaMethods)),
						'MA' = array(rep(metaMethods,each=numpoolmeths*number.thresholds)),
						'source' = rep(names_one,each=number.thresholds))

# Calculate and sort by AUC values (area under curve): random, fixed and ALE
LFPR <- split(MEANFPR, col(MEANFPR))
LTPR <- split(MEANTPR, col(MEANTPR))
AUC_VALUES <- mapply(LFPR, FUN = comp_auc, TP = LTPR)
	names(AUC_VALUES) <- names_one
	TFPR.DATA$MA <- factor(TFPR.DATA$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
	TFPR.DATA$pooling <- factor(TFPR.DATA$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = ArtPOOLlabels)

# What is the point in number.thresholds most close to alpha = 0.05
thresholds <- seq(0,1,length.out = number.thresholds)
	IDthreshold <- which(abs(c(thresholds - 0.05)) == min(abs(c(thresholds - 0.05))))

lineStarts <- lineEnds <- c()
for(i in 1:length(names(FPR))){
	line <- TFPR.DATA[TFPR.DATA$source == names_one[i],'FP'][IDthreshold]
		lineStarts <- c(lineStarts, line)
	line_end <- TFPR.DATA[TFPR.DATA$source == names_one[i],'TP'][IDthreshold]
		lineEnds <- c(lineEnds, line_end)
}
	TFPR.DATA$LineStart <- rep(lineStarts, each = number.thresholds)
	TFPR.DATA$LineEnd <- rep(lineEnds, each = number.thresholds)

# All combinations of pooling and MA
combinations <- data.frame(expand.grid(
					factor(poolmeth,levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = c('OLS', 'Fixed Effects', 'Mixed Effects')),
					factor(metaMethods, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))))
colnames(combinations) <- c('pooling', 'MA')
AUC_DAT <- data.frame(combinations,
			AUC = format(round(unlist(AUC_VALUES),4), nsmall = 4))


##
###############
### Scaled partial AUC
###############
##

## We need to reshape FPR and TPR to long format
# Start with going from list to data frame
FPR_long <- do.call(rbind.data.frame, FPR) %>%
		# Bind rownames to the data frame: we will use these for names of MA and pooling method
		bind_cols(.,data.frame(source = rownames(.))) %>%
			# Now we only take this column and extract the name of the pooling method
			mutate(., pooling = regmatches(source, regexpr(paste(poolmeth, collapse = "|"),source))) %>%
			# Same for MA
			mutate(., MA = regmatches(source, regexpr(paste(metaMethods, collapse = "|"),source))) %>%
			# Drop this column
			select(., -source) %>%
			# Reshape to long format
			gather(., key = Run, value = FPR, -pooling, -MA)

## Same for TPR
# Start with going from list to data frame
TPR_long <- do.call(rbind.data.frame, TPR) %>%
		# Bind rownames to the data frame: we will use these for names of MA and pooling method
		bind_cols(.,data.frame(source = rownames(.))) %>%
			# Now we only take this column and extract the name of the pooling method
			mutate(., pooling = regmatches(source, regexpr(paste(poolmeth, collapse = "|"),source))) %>%
			# Same for MA
			mutate(., MA = regmatches(source, regexpr(paste(metaMethods, collapse = "|"),source))) %>%
			# Drop this column
			select(., -source) %>%
			# Reshape to long format
			gather(., key = Run, value = TPR, -pooling, -MA)

# Combine TPR and FPR
RATES_long <- bind_cols(FPR_long, data.frame(TPR = TPR_long$TPR, NSTUD = NSTUD)) %>% tbl_df()

# Gather AUC values in long format
StPartAUC <- RATES_long %>%  group_by(pooling, MA, Run) %>% slice(1:11) %>%
				do(PARTauc = comp_auc(FP=.$FPR, TP = .$TPR)) %>%
				mutate(., PARTauc = unlist(PARTauc)) %>% tbl_df() %>%
				mutate(., sPARTauc = 1/2*(1 + (PARTauc - 0.1**2/2) / (0.1 - 0.1**2/2))) %>%
				ungroup() %>% group_by(pooling, MA) %>%
				summarise(AvSPartauc = mean(sPARTauc)) %>%
				mutate(AvSPartauc = round(AvSPartauc,4))

# Label the factors
StPartAUC$pooling <- factor(StPartAUC$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = c('OLS', 'Fixed Effects', 'Mixed Effects'))
StPartAUC$MA <- factor(StPartAUC$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))


##
###############
### Plotting
###############
##

# Number of ROC figures in paper
ROCFigNumber <- c(3,5,7)

quartz(height = 6.5, width = 3.5,
	type = 'png',file = paste(LocFileSave, '/figure',ROCFigNumber[WD],'_ROC_K',NSTUD,'.png', sep = ''), dpi = 600)
	ggplot(TFPR.DATA, aes(FP,TP)) +
		geom_point(colour = "black", size = 0.05) +
		scale_colour_manual(guide = FALSE) +
		facet_wrap(MA ~ pooling, dir = "v") +
    geom_ribbon(aes(ymin = TP-SDTP, ymax = TP+SDTP), alpha=0.2) +
		geom_segment(aes(x = LineStart, y = 0, xend = LineStart, yend = LineEnd), colour = "black", linetype = 2 , size = 0.2) +
		geom_segment(aes(x = 0, y = LineEnd, xend = LineStart, yend = LineEnd), colour = "black", linetype = 2, size = 0.2) +
		geom_text(aes(x = 0.70, y = 0.19, label=AUC), data=AUC_DAT, size = 3.0) +
		scale_x_continuous(name = paste('False Positive Rate (K = ', NSTUD,')', sep = '')) +
		scale_y_continuous(name = paste('True Positive Rate (K = ', NSTUD,')', sep = '')) +
		theme_bw(base_size = 9, base_family = "Helvetica") +
		theme(legend.position='bottom',
					legend.text = element_text(size=2.5),
					legend.title = element_text(size=2.5),
					panel.grid.major = element_line(size = 0.2, color = "grey"),
					axis.line = element_line(size = 0.4, color = "black"),
					axis.title.y = element_text(size=9.2),
					axis.title.x = element_text(size=9.2),
					axis.text.x = element_text(angle = 45, vjust = 0.25,hjust = 0.5),
					strip.background = element_rect(colour = "white", fill = "white"))
dev.off()


# Number of partial ROC figures in paper
partROCFigNumber <- c(4,6,8)

# Partial ROC curve: from alpha 0 - 0.1
quartz(height = 6.5, width = 3.5,
	type = 'png',file = paste(LocFileSave, '/figure',partROCFigNumber[WD],'_PartROC_K',NSTUD,'.png', sep = ''), dpi = 600)
TFPR.DATA  %>%
	group_by(source) %>% slice(.,1:11) %>%
		ggplot(., aes(FP,TP)) +
		geom_point(size = 0.3) +
		geom_line(size = 0.45) +
		facet_wrap(MA ~ pooling, dir = "v") +
    geom_ribbon(aes(ymin = TP-SDTP, ymax = TP+SDTP), alpha=0.2) +
		geom_segment(aes(x = LineStart, y = 0, xend = LineStart, yend = LineEnd), linetype = 2 , size = 0.25) +
		geom_segment(aes(x = 0, y = LineEnd, xend = LineStart, yend = LineEnd), linetype = 2, size = 0.25) +
		scale_x_continuous(name = paste('False Positive Rate (K = ', NSTUD,')', sep = ''), breaks = seq(0,0.1,by = 0.02)) +
		scale_y_continuous(name = paste('True Positive Rate (K = ', NSTUD,')', sep = ''), limits = c(0,0.75)) +
		geom_text(aes(x = 0.024, y = 0.72, label=AvSPartauc), data=StPartAUC, size = 3.05) +
		theme_bw(base_size = 9, base_family = "Helvetica") +
		theme(legend.position='bottom',
					legend.text = element_text(size=9),
					legend.title = element_text(size=9),
					panel.grid.major = element_line(size = 0.2, color = "grey"),
					axis.line = element_line(size = 0.4, color = "black"),
					axis.title.y = element_text(size=9.2),
					axis.title.x = element_text(size=9.2),
					axis.text.x = element_text(angle = 45, vjust = 0.25,hjust = 0.5),
					strip.background = element_rect(colour = "white", fill = "white"))
dev.off()









################################################################################
################################################################################
################################### APPENDIX ###################################
################################################################################
################################################################################

##
###############
### Plotting AUC and partial AUC for FDR 0.05
###############
##

# For efficiency, I will re-use objects. So better remove them first.
rm(FPR.DATA,TPR.DATA,SDTPR.DATA,TFPR.DATA,LFPR,
	LTPR,AUC_VALUES,thresholds,IDthreshold,
	lineStarts,lineEnds,combinations,AUC_DAT)

# Mean and sd
MEANFPR.tmp_FDR05 <- lapply(FPR_FDR05,rowMeans)
	MEANFPR_FDR05 <- t(do.call(rbind, MEANFPR.tmp_FDR05))
MEANTPR.tmp_FDR05 <- lapply(TPR_FDR05,rowMeans)
SDTPR.tmp_FDR05 <- lapply(TPR_FDR05, function(x) apply(x, 1, sd))
	MEANTPR_FDR05 <- t(do.call(rbind, MEANTPR.tmp_FDR05))
	SDTPR_FDR05 <- t(do.call(rbind, SDTPR.tmp_FDR05))

# Put them in data frames
FPR.DATA <- data.frame('value' = matrix(MEANFPR_FDR05,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

TPR.DATA <- data.frame('value' = matrix(MEANTPR_FDR05,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

SDTPR.DATA<- data.frame('value' = matrix(SDTPR_FDR05,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

TFPR.DATA <- data.frame('FP' = matrix(MEANFPR_FDR05,ncol=1),
						'TP' = matrix(MEANTPR_FDR05,ncol=1),
						'SDTP' = matrix(SDTPR_FDR05,ncol=1),
						'pooling' = array(rep(rep(poolmeth,each=number.thresholds),numMetaMethods)),
						'MA' = array(rep(metaMethods,each=numpoolmeths*number.thresholds)),
						'source' = rep(names_one,each=number.thresholds))

# Calculate and sort by AUC values (area under curve): random, fixed and ALE
LFPR <- split(MEANFPR_FDR05, col(MEANFPR_FDR05))
LTPR <- split(MEANTPR_FDR05, col(MEANTPR_FDR05))
AUC_VALUES <- mapply(LFPR, FUN = comp_auc, TP = LTPR)
	names(AUC_VALUES) <- names_one
	TFPR.DATA$MA <- factor(TFPR.DATA$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
	TFPR.DATA$pooling <- factor(TFPR.DATA$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = ArtPOOLlabels)

# What is the point in number.thresholds most close to alpha = 0.05
thresholds <- seq(0,1,length.out = number.thresholds)
	IDthreshold <- which(abs(c(thresholds - 0.05)) == min(abs(c(thresholds - 0.05))))

lineStarts <- lineEnds <- c()
for(i in 1:length(names(FPR))){
	line <- TFPR.DATA[TFPR.DATA$source == names_one[i],'FP'][IDthreshold]
		lineStarts <- c(lineStarts, line)
	line_end <- TFPR.DATA[TFPR.DATA$source == names_one[i],'TP'][IDthreshold]
		lineEnds <- c(lineEnds, line_end)
}
	TFPR.DATA$LineStart <- rep(lineStarts, each = number.thresholds)
	TFPR.DATA$LineEnd <- rep(lineEnds, each = number.thresholds)

# All combinations of pooling and MA
combinations <- data.frame(expand.grid(
					factor(poolmeth,levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = c('OLS', 'Fixed Effects', 'Mixed Effects')),
					factor(metaMethods, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))))
colnames(combinations) <- c('pooling', 'MA')
AUC_DAT <- data.frame(combinations,
			AUC = format(round(unlist(AUC_VALUES),4), nsmall = 4))

# Mean AUC for the MA's
meanAUCperMA <- aggregate(as.numeric(as.character(AUC)) ~ MA, FUN = mean, data = AUC_DAT)
meanAUCperMA


ggplot(TFPR.DATA, aes(FP,TP)) +
		geom_point(colour = "black", size = 0.7) +
		scale_colour_manual(values = coloursROC, guide = FALSE) +
		facet_wrap(MA ~ pooling, dir = "v") +
    geom_ribbon(aes(ymin = TP-SDTP, ymax = TP+SDTP), alpha=0.2) +
		geom_segment(aes(x = LineStart, y = 0, xend = LineStart, yend = LineEnd), colour = "black", linetype = 2 , size = 0.5) +
		geom_segment(aes(x = 0, y = LineEnd, xend = LineStart, yend = LineEnd), colour = "black", linetype = 2, size = 0.5) +
		geom_text(aes(x = 0.82, y = 0.30, label=AUC), data=AUC_DAT, size = 4) +
		scale_x_continuous(name = paste('False Positive Rate (K = ', NSTUD,')', sep = '')) +
		scale_y_continuous(name = paste('True Positive Rate (K = ', NSTUD,')', sep = '')) +
		theme_bw(base_size = 10, base_family = "Helvetica") +
		ggtitle('FDR 0.05') +
		theme(legend.position='bottom',
					legend.text = element_text(size=9),
					legend.title = element_text(size=9),
					panel.grid.major = element_line(size = 0.2, color = "grey"),
					axis.line = element_line(size = 0.8, color = "black"),
					axis.title.y = element_text(size=9.2),
					axis.title.x = element_text(size=9.2),
					strip.background = element_rect(colour = "white", fill = "white"))

ggsave(file = paste(LocFileSave, '/sup_ROC_FDR05_K',NSTUD,'.pdf', sep = ''), device = "pdf",
			width = 4.8, height = 7, units = 'in', dpi = 300)

##
###############
### Plotting AUC for uncorrected 0.001
###############
##

# For efficiency, I will re-use objects. So better remove them first.
rm(FPR.DATA,TPR.DATA,SDTPR.DATA,TFPR.DATA,LFPR,
	LTPR,AUC_VALUES,thresholds,IDthreshold,
	lineStarts,lineEnds,combinations,AUC_DAT,meanAUCperMA)

# Mean and sd
MEANFPR.tmp_Un001 <- lapply(FPR_Un001,rowMeans)
	MEANFPR_Un001 <- t(do.call(rbind, MEANFPR.tmp_Un001))
MEANTPR.tmp_Un001 <- lapply(TPR_Un001,rowMeans)
SDTPR.tmp_Un001 <- lapply(TPR_Un001, function(x) apply(x, 1, sd))
	MEANTPR_Un001 <- t(do.call(rbind, MEANTPR.tmp_Un001))
	SDTPR_Un001 <- t(do.call(rbind, SDTPR.tmp_Un001))

# Put them in data frames
FPR.DATA <- data.frame('value' = matrix(MEANFPR_Un001,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

TPR.DATA <- data.frame('value' = matrix(MEANTPR_Un001,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

SDTPR.DATA<- data.frame('value' = matrix(SDTPR_Un001,ncol=1),
						'source' = rep(names_one,each=number.thresholds))

TFPR.DATA <- data.frame('FP' = matrix(MEANFPR_Un001,ncol=1),
						'TP' = matrix(MEANTPR_Un001,ncol=1),
						'SDTP' = matrix(SDTPR_Un001,ncol=1),
						'pooling' = array(rep(rep(poolmeth,each=number.thresholds),numMetaMethods)),
						'MA' = array(rep(metaMethods,each=numpoolmeths*number.thresholds)),
						'source' = rep(names_one,each=number.thresholds))

# Calculate and sort by AUC values (area under curve): random, fixed and ALE
LFPR <- split(MEANFPR_Un001, col(MEANFPR_Un001))
LTPR <- split(MEANTPR_Un001, col(MEANTPR_Un001))
AUC_VALUES <- mapply(LFPR, FUN = comp_auc, TP = LTPR)
	names(AUC_VALUES) <- names_one
	TFPR.DATA$MA <- factor(TFPR.DATA$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
	TFPR.DATA$pooling <- factor(TFPR.DATA$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = ArtPOOLlabels)

# What is the point in number.thresholds most close to alpha = 0.05
thresholds <- seq(0,1,length.out = number.thresholds)
	IDthreshold <- which(abs(c(thresholds - 0.05)) == min(abs(c(thresholds - 0.05))))

lineStarts <- lineEnds <- c()
for(i in 1:length(names(FPR))){
	line <- TFPR.DATA[TFPR.DATA$source == names_one[i],'FP'][IDthreshold]
		lineStarts <- c(lineStarts, line)
	line_end <- TFPR.DATA[TFPR.DATA$source == names_one[i],'TP'][IDthreshold]
		lineEnds <- c(lineEnds, line_end)
}
	TFPR.DATA$LineStart <- rep(lineStarts, each = number.thresholds)
	TFPR.DATA$LineEnd <- rep(lineEnds, each = number.thresholds)

# All combinations of pooling and MA
combinations <- data.frame(expand.grid(
					factor(poolmeth,levels = c('OLS', 'FixedEffect', 'MixedEffect'), labels = c('OLS', 'Fixed Effects', 'Mixed Effects')),
					factor(metaMethods, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))))
colnames(combinations) <- c('pooling', 'MA')
AUC_DAT <- data.frame(combinations,
			AUC = format(round(unlist(AUC_VALUES),4), nsmall = 4))

# Mean AUC for the MA's
meanAUCperMA <- aggregate(as.numeric(as.character(AUC)) ~ MA, FUN = mean, data = AUC_DAT)
meanAUCperMA


ggplot(TFPR.DATA, aes(FP,TP)) +
		geom_point(colour = "black", size = 0.7) +
		scale_colour_manual(values = coloursROC, guide = FALSE) +
		facet_wrap(MA ~ pooling, dir = "v") +
    geom_ribbon(aes(ymin = TP-SDTP, ymax = TP+SDTP), alpha=0.2) +
		geom_segment(aes(x = LineStart, y = 0, xend = LineStart, yend = LineEnd), colour = "black", linetype = 2 , size = 0.5) +
		geom_segment(aes(x = 0, y = LineEnd, xend = LineStart, yend = LineEnd), colour = "black", linetype = 2, size = 0.5) +
		geom_text(aes(x = 0.82, y = 0.30, label=AUC), data=AUC_DAT, size = 4) +
		scale_x_continuous(name = paste('False Positive Rate (K = ', NSTUD,')', sep = '')) +
		scale_y_continuous(name = paste('True Positive Rate (K = ', NSTUD,')', sep = '')) +
		theme_bw(base_size = 10, base_family = "Helvetica") +
		ggtitle('Uncorrected 0.001') +
		theme(legend.position='bottom',
					legend.text = element_text(size=9),
					legend.title = element_text(size=9),
					panel.grid.major = element_line(size = 0.2, color = "grey"),
					axis.line = element_line(size = 0.8, color = "black"),
					axis.title.y = element_text(size=9.2),
					axis.title.x = element_text(size=9.2),
					strip.background = element_rect(colour = "white", fill = "white"))

ggsave(file = paste(LocFileSave, '/sup_ROC_Un001_K',NSTUD,'.pdf', sep = ''), device = "pdf",
			width = 4.8, height = 7, units = 'in', dpi = 300)




