####################
#### TITLE:     Calculate FPR/TPR when comparing MA's with a reference from the test condition.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/Analyses
#### First Modified: 21/11/2017
#### Notes:
#################

# Notes:
# This file is used to gather the raw data and process them to intermediate.
# These are then saved in this directory.

# - # We process:
# --- # data frame with TPR/FPR
# --- # data frame with overlap values 

# Figures for appendix are created here, as otherwise I need to save too much intermediate files...

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


# Results, in list according to K (number of studies in MA)
WDs <- list(
  '10' = "/<<LOCATION_OF_RESULTS>>/K10",
  '12' = "/<<LOCATION_OF_RESULTS>>/K12",
  '14' = "/<<LOCATION_OF_RESULTS>>/K14",
  '16' = "/<<LOCATION_OF_RESULTS>>/K16",
  '18' = "/<<LOCATION_OF_RESULTS>>/K18",
  '20' = "/<<LOCATION_OF_RESULTS>>/K20",
  '30' = "/<<LOCATION_OF_RESULTS>>/K30",
  '35' = "/<<LOCATION_OF_RESULTS>>/K35"
)

# Source the directories from Directories.R. This file contains the paths to the raw data.
# It is added to the .gitignore file.
if(!grepl(pattern = 'IntermediateResults', x = getwd())){
  stop('Working directory not at IntermediateResults')
}
source('../Directories.R')

# Location to save figures for appendix
LocFileSave <- '../Figures/Appendix'

##
###############
### Functions/libraries
###############
##

# Libraries
library(ggplot2)
library(ggthemes)
library(oro.nifti)
library(dplyr)
library(tidyr)

# Function to calculate ROC values based on P-map
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

# Function that calculates auc on columns of FPR and TPR
returnAUC <- function(FPR, TPR){
  t(mapply(comp_auc, split(FPR, col(FPR)), 
           TP = split(TPR, col(TPR))))
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


##
###############
### Global values
###############
##

# Choose your WD
WD <- 1

# Data frame with:
# -- Working directories (WDs) to choose from
# -- number of runs/folds (see the SamplingCVContK.R file for the numbers)
# -- number of studies in the MA (K) and 
# -- number of subjects in the reference image.
DesignInfo <- data.frame(WDs = 1:8,
                FOLDS = c(7, 5, 5, 4, 3, 3, 2, 2),
                K = c(10, 12, 14, 16, 18, 20, 30, 35),
                NSUBREF = c(200, 240, 280, 320, 360, 400, 600, 700))
print(DesignInfo)
# Now select NRUNS, NSTUD
NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
  unlist() %>% as.numeric()
NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
  unlist() %>% as.numeric()

# Dimension of the original data
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

################################################################################
################################# TPR and FPR ##################################
################################################################################

##
###############
### TPR/FPR values
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

##
###############
### Calculate TPR/FPR
###############
##

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

# Add the number of studies in the MA to the data frame
TFPR.DATA$K <- NSTUD

# Calculate and sort by AUC values (area under curve): random, fixed and ALE
LFPR <- split(MEANFPR, col(MEANFPR))
LTPR <- split(MEANTPR, col(MEANTPR))
AUC_VALUES <- mapply(LFPR, FUN = comp_auc, TP = LTPR)
  names(AUC_VALUES) <- names_one

# Calculate AUC on individual FP and TP of each run
AUC_long <- data.frame() %>% tbl_df()
for(i in 1:length(FPR)){
  unlFPR <- FPR[[i]]
  unlTPR <- TPR[[i]]
  AUC_long <- data.frame('AUC' = as.numeric(returnAUC(FPR = unlFPR, TPR = unlTPR))) %>%
    mutate(Source = names(FPR)[i],
           pooling = names_two[i, 'Pooling'],
           MA = names_two[i, 'MetaAnalysis'],
           Run = 1:NRUNS,
           K = NSTUD) %>% tbl_df() %>%
    bind_rows(AUC_long, .)
}


##
###############
### Write R objects to folder
###############
##

# # False positive rates
# saveRDS(FPR, file = paste('FPR/FPR_K_', NSTUD,'.rds', sep = ''))
# 
# # True positive rates
# saveRDS(TPR, file = paste('TPR/TPR_K_', NSTUD,'.rds', sep = ''))

# Data frame with false and true positive rates
saveRDS(TFPR.DATA, file = paste('TFPR/TFPR_K_', NSTUD,'.rds', sep = ''))

# Data frames with summary of AUC and AUC in long format
saveRDS(AUC_VALUES, file = paste('AUC_avgFPTP/AUC_VALUES_K_', NSTUD,'.rds', sep = ''))
saveRDS(AUC_long, file = paste('AUC_long/AUC_long_K_', NSTUD,'.rds', sep = ''))


################################################################################
################################################################################
################################### APPENDIX ###################################
################################################################################
################################################################################

APPENDIX <- FALSE

# Run if necessary
if(APPENDIX){

# Labels of appendix figures in paper:
# -- 1 - 3: partial ROC for K = 10, 20 and 35
# -- 4 - 6: FDR = 0.05 for K = 10, 20 and 35
# -- 7 - 9: Uncorrected = 0.001 for K = 10, 20 and 35

AppLabels <- c(rep('A2', 3), rep('A3', 3), rep('A4', 3))



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
StPartAUC <- RATES_long %>%  group_by(pooling, MA, Run) %>% dplyr::slice(1:11) %>%
  do(PARTauc = comp_auc(FP=.$FPR, TP = .$TPR)) %>%
  mutate(., PARTauc = unlist(PARTauc)) %>% tbl_df() %>%
  mutate(., sPARTauc = 1/2*(1 + (PARTauc - 0.1**2/2) / (0.1 - 0.1**2/2))) %>%
  ungroup() %>% group_by(pooling, MA) %>%
  summarise(AvSPartauc = mean(sPARTauc)) %>%
  mutate(AvSPartauc = round(AvSPartauc,4))

# Label the factors
StPartAUC$pooling <- factor(StPartAUC$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), 
                            labels = c('OLS', 'Fixed Effects', 'Mixed Effects'))
StPartAUC$MA <- factor(StPartAUC$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), 
                       labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
TFPR.DATA$MA <- factor(TFPR.DATA$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), 
                       labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
TFPR.DATA$pooling <- factor(TFPR.DATA$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), 
                            labels = ArtPOOLlabels)

# What is the point in number.thresholds most close to alpha = 0.05
thresholds <- seq(0,1,length.out = number.thresholds)
IDthreshold <- which(abs(c(thresholds - 0.05)) == min(abs(c(thresholds - 0.05))))

# Calculate the start and end of the line drops: most closest to IDthreshold
lineStarts <- lineEnds <- c()
for(i in 1:length(names_one)){
  line <- TFPR.DATA[TFPR.DATA$source == names_one[i],'FP'][IDthreshold]
  lineStarts <- c(lineStarts, line)
  line_end <- TFPR.DATA[TFPR.DATA$source == names_one[i],'TP'][IDthreshold]
  lineEnds <- c(lineEnds, line_end)
}

# Add to data frame
TFPR.DATA$LineStart <- rep(lineStarts, each = number.thresholds)
TFPR.DATA$LineEnd <- rep(lineEnds, each = number.thresholds)

# Partial ROC curve: from alpha 0 - 0.1
quartz(height = 6.5, width = 3.5,
       type = 'png',file = paste(LocFileSave, '/sup_',AppLabels[WD],'_PartROC_K',NSTUD,'.png', sep = ''), dpi = 600)
TFPR.DATA  %>%
  group_by(source) %>% dplyr::slice(.,1:11) %>%
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

ggsave(file = paste(LocFileSave, '/sup_', AppLabels[WD + 3], '_ROC_FDR05_K',NSTUD,'.pdf', sep = ''), device = "pdf",
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

ggsave(file = paste(LocFileSave, '/sup_', AppLabels[WD+6], '_ROC_Un001_K',NSTUD,'.pdf', sep = ''), device = "pdf",
       width = 4.8, height = 7, units = 'in', dpi = 300)

# End of appendix if statement
}


################################################################################
################################### OVERLAP ####################################
################################################################################

##
###############
### Overlap preparation
###############
##

# Mask for thresholding the PMaps of fixed and random effects MA: note that we used one universal mask throughout the entire study.
MASK <- readNIfTI(paste(WDs[[WD]],'/Run_1/GroupAnalysis/mask.nii.gz',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
# Mask for ALE
MASKALE <- MNI152

# Array of meta-analysis methods: different than TFPR
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

##
###############
### Write R object to folder
###############
##

saveRDS(OverlapP, file = paste('Overlap/Overlap_K_', NSTUD,'.rds', sep = ''))




################################################################################
################################################################################
################################### APPENDIX ###################################
################################################################################
################################################################################



APPENDIX <- FALSE

# Run if necessary
if(APPENDIX){
  # Label for figure  
  AppLabels <- paste('A', 6, sep = '')
  
  ##
  ###############
  ### Calculate overlap for appendix: uncorrected 0.001
  ###############
  ##
  
  # Array of meta-analysis methods for overlap in appendix
  metaMethodsApp <- c('FixedUnApp', 'RanUnApp', 'ALEUnApp')
  numMetaMethodsApp <- length(metaMethodsApp)
  
  # number of columns in the dataset for the appendix
  numColsApp <- numMetaMethodsApp * numpoolmeths
  
  # List of all the names for loading in the data
  METHODSapp <- list(
    array(rep(metaMethodsApp,each=numpoolmeths)),
    array(rep(poolmeth,numMetaMethodsApp)))
  names(METHODSapp) <- c('MetaAnalysis', 'Pooling')
  
  # Original names to load in the data files
  METHODS <- list(
    array(rep(metaMethods,each=numpoolmeths)),
    array(rep(poolmeth,numMetaMethods)))
  names(METHODS) <- c('MetaAnalysis', 'Pooling')
  
  # Labels for meta-analyses
  ArtLABMA_app <- c('Fixed Effects MA', 'Random Effects MA','ALE Uncorrected')
  
  # Empty data frame
  Overlap_KApp <- data.frame()
  
  # For loop over the WDs 
  for(w in 1:length(WDs)){
    print(paste0(' At WD: ', w))
    # Corresponding WD
    WD <- w
    
    # Now select NRUNS, NSTUD
    NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
      unlist() %>% as.numeric()
    NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
      unlist() %>% as.numeric()
    
    # The pairwise comparissons
    combRuns <- c(1:NRUNS)
    PAIRS <- t(combn(combRuns,2))
    NPAIRS <- dim(PAIRS)[1]
    
    # Overlap in each PAIR
    OverlapApp <- list(
      array(0,dim=c(numpoolmeths,NPAIRS),dimnames=list(poolmeth,paste('P',c(1:NPAIRS),sep=''))),
      array(0,dim=c(numpoolmeths,NPAIRS),dimnames=list(poolmeth,paste('P',c(1:NPAIRS),sep=''))),
      array(0,dim=c(numpoolmeths,NPAIRS),dimnames=list(poolmeth,paste('P',c(1:NPAIRS),sep='')))
    )
    names(OverlapApp) <- metaMethodsApp
    
    # For loop over the PAIRS
    for(p in 1:NPAIRS){
      # At run:
      print(paste('At pair ',p,sep=''))
      # Now we loop over the pooling and meta-analysis methods
      # -- conveniently: we use the METHODS object from main text. Here, ALECluster is the last object (position 10-12),
      # so by looping only to first 9 objects, we can re-use METHODS and skip ALECluster.
      for(j in 1:numColsApp){
        if(grepl('ALE', METHODSapp[['MetaAnalysis']][j])){
          # Assign the MASKALE to the variable OVERLAPMASK: used when calculating overlap
          assign('OVERLAPMASK', MASKALE)
          # Load in the data: ALE uncorrected using matlab code of Eickhoff, then cluster corrected
          if(METHODSapp[['MetaAnalysis']][j] == "ALEUnApp"){
            # Group A: Z-values then going to P-values
            GroupA.Z <- readNIfTI(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/ALE/',METHODS[[2]][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
            GroupA.tmp <- 1-pnorm(GroupA.Z, mean = 0, sd = 1)
            
            # Group B
            GroupB.Z <- readNIfTI(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/ALE/',METHODSapp[[2]][j],'/ALE/ALEvolumesZ/OLS.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
            GroupB.tmp <- 1-pnorm(GroupB.Z, mean = 0, sd = 1)
            
            # Thresholding. Maps are already thresholded, hence NA to argument.
            GroupA <- ThreshPVal(GroupA.tmp, 0.001, mask = MASKALE)
            GroupB <- ThreshPVal(GroupB.tmp, 0.001, mask = MASKALE)
            rm(GroupA.tmp,GroupB.tmp,GroupA.Z,GroupB.Z)
          }
        }else{
          # Assign the MASK to the variable OVERLAPMASK
          assign('OVERLAPMASK', MASK)
          # Load in PVal of group A: fixed and random effects MA
          load(paste(WDs[[WD]],'/Run_',PAIRS[p,1],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/PVal',sep=''))
          
          # Analyses for appendix: only test at P < 0.001
          GroupA <- ThreshPVal(PVal, threshold = 0.001, mask = MASK)
          rm(PVal)
  
          # Load in PVal of group B
          load(paste(WDs[[WD]],'/Run_',PAIRS[p,2],'/MetaAnalyses/',METHODS[[1]][j],'/',METHODS[[2]][j],'/PVal',sep=''))
          GroupB <- ThreshPVal(PVal, threshold = 0.001, mask = MASK)
          rm(PVal)
        }
        # Now that we have the thresholded map for group A and B, let us calculate the pairwise overlap
        method <- METHODSapp[[1]][j]
        pooling <- METHODSapp[[2]][j]
        OverlapApp[[method]][pooling,p] <- OVERLAP(GroupA, GroupB, mask = OVERLAPMASK)
      }
      # Remove masks
      rm(OVERLAPMASK)
    }
  
    # Now go from list to data frame
    OverlapDApp <- data.frame('Overlap' = unlist(OverlapApp),
                           'pooling' = rep(poolmeth, NPAIRS * numMetaMethodsApp),
                           'MA' = rep(names(OverlapApp), each = NPAIRS * numpoolmeths)) %>%
      `rownames<-`(1:c(NPAIRS * numMetaMethodsApp * numpoolmeths)) %>%
      mutate(K = NSTUD, Pair = rep(1:NPAIRS, numMetaMethodsApp * numpoolmeths))
    
    # Add to data frame
    Overlap_KApp <- bind_rows(Overlap_KApp, OverlapDApp)
  }
  
  # Add labeling for plotting to factors
  Overlap_KApp$pooling <- factor(Overlap_KApp$pooling, levels = sort(poolmeth), labels = sort(ArtPOOLlabels))
  Overlap_KApp$MA <- factor(Overlap_KApp$MA, levels = sort(metaMethodsApp), labels = sort(ArtLABMA_app))
  
  # Plot in figure 
  Overlap_KApp %>% group_by(pooling, MA, K) %>% 
    summarise(AvgOv = mean(Overlap)) %>% 
    ggplot(aes(x = K, y = AvgOv)) +
    geom_line(aes(colour = pooling), size = 0.9) +
    geom_point(aes(colour = pooling), size = 0.9) + 
    facet_grid(~ MA) +
    scale_x_continuous(name = 'Number of studies (K)') +
    scale_y_continuous(name = 'Average overlap of activation') +
    scale_color_brewer('Group level model', type = 'qual', palette = 2) +
    ggtitle('Average overlap for each MA thresholded at:', 
            subtitle =  'P (uncorrected) < 0.001') +
    scale_alpha_manual('', values = c(1,0.5), guide = 'none') + 
    theme_bw(base_size = 9, base_family = "Helvetica") +
    theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
          axis.line = element_line(size = 0.4, color = "black"),
          axis.title.y = element_text(size=11),
          axis.text.y = element_text(size=10),
          axis.title.x = element_text(size=11),
          axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
          strip.background = element_rect(colour = "white", fill = "white"))

  # Save figure
  ggsave(file = paste(LocFileSave, '/sup_', AppLabels[1], '_AvgOverlap_Un001.pdf', sep = ''), device = "pdf",
         height = 5, width = 8.2, units = 'in', dpi = 300)
}
