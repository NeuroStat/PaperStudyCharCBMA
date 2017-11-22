####################
#### TITLE:     Calculate ROC curves when comparing MA's with a reference from the test condition.
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

# Check if we are in the correct working directory
if(!grepl(pattern = 'Analyses', x = getwd())){
  stop('Working directory not at IntermediateResults')
}

# MNI152 template: will be used for masking in ALE
#MNI152 <- readNIfTI('<<LOCATION_TO_MNI152_TEMPLATE>>/MNI152.nii')[,,]
#	IDMNI152 <- MNI152 == 0

# Location to save plots for the paper
LocFileSave <- '/Figures/'

# Results, in list according to K (number of studies in MA)
WDs <- list(
	'10' = "IntermediateResults",
	'20' = "/<<LOCATION_OF_RESULTS>>/K20",
	'35' = "/<<LOCATION_OF_RESULTS>>/K35"
	)

# Libraries
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)

##
###############
### Preparation
###############
##

# Choose your WD
WD <- 1

# Data frame with:
# -- Possible working directories (WDs)
# -- number of runs/folds
# -- number of studies in the MA (K) and 
# -- number of subjects in the reference image.
DesignInfo <- data.frame(WDs = 1:3,
                         FOLDS = c(7,3,2),
                         K = c(10,20,35),
                         NSUBREF = c(200,400,700))
print(DesignInfo)

# Vector of pooling methods
poolmeth <- c(
  'FixedEffect',
  'OLS',
  'MixedEffect'
)
numpoolmeths <- length(poolmeth)

# Array of meta-analysis methods
metaMethods <- c('FixedEffUn', 'RanEffUn','ALE')
numMetaMethods <- length(metaMethods)

# List of all the pooling methods and meta-analyses when loading data
METHODS <- list(
  array(rep(metaMethods,each=numpoolmeths)),
  array(rep(poolmeth,numMetaMethods)))
names(METHODS) <- c('MetaAnalysis', 'Pooling')

# Names of objects: CARE!!
# The order = each MA first, then the pooling methods. 
# This corresponds with METHODS object.
names_two <- unique(expand.grid(rev(METHODS)))
names_one <- paste(names_two[,1],names_two[,2], sep=":")

# Labels to be used in paper for group level models
ArtPOOLlabels <- c('OLS', 'Fixed Effects', 'Mixed Effects')

# Labels of meta-analyses
ArtLABMA <- c('ALE', 'Fixed Effects MA','Random Effects MA')

# number of thresholds considered in the ROC
number.thresholds <- 100


##
###############
### Load values for FP, TP and AUC + calculate line drops
###############
##

# Empty TFPR and AUC data frames
TFPR.DATA_K <- AUC_K <- data.frame() %>% as_tibble()

for(f in 1:length(WDs)){
  # WD gets f in loop 
  WD <- f
  
  # Select NRUNS, NSTUD
  NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
    unlist() %>% as.numeric()
  NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
    unlist() %>% as.numeric()
  
  # Loop over the number of data frames to load in
  TFPR.DATA <- readRDS(paste('IntermediateResults/TFPR_K_', NSTUD,'.rds', sep = ''))
  AUC_VALUES <- readRDS(paste('IntermediateResults/AUC_VALUES_K_', NSTUD,'.rds', sep = ''))
  
  # First re-label the factors
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
  
  # Add to data frame, then convert to tibble
  TFPR.DATA$LineStart <- rep(lineStarts, each = number.thresholds)
  TFPR.DATA$LineEnd <- rep(lineEnds, each = number.thresholds)
  TFPR.DATA <- as_tibble(TFPR.DATA)
  
  # Create all combinations of pooling and MA for labels of AUC
  combinations <- data.frame(expand.grid(
    factor(poolmeth,levels = c('OLS', 'FixedEffect', 'MixedEffect'), 
           labels = c('OLS', 'Fixed Effects', 'Mixed Effects')),
    factor(metaMethods, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), 
           labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))))
  colnames(combinations) <- c('pooling', 'MA')
  AUC_DAT <- data.frame(combinations,
                AUC = AUC_VALUES %>% unlist() %>% 
                    round(.,4) %>% format(., nsmall = 4) %>% as.numeric())
  
  # Add K to data frame and convert to tibble
  AUC_DAT <- AUC_DAT %>% mutate(K = NSTUD) %>% as_tibble()

  # Bind TFPR and AUC to data frames
  TFPR.DATA_K <- bind_rows(TFPR.DATA_K, TFPR.DATA)
  AUC_K <- bind_rows(AUC_K, AUC_DAT)
}

# 
# # Loop over the number of data frames to load in
# TFPR.DATA <- readRDS(paste('IntermediateResults/TFPR_K_', NSTUD,'.rds', sep = ''))
# AUC_VALUES <- readRDS(paste('IntermediateResults/AUC_VALUES_K_', NSTUD,'.rds', sep = ''))
# 
# 
# 
# ##
# ###############
# ### Calculate line drops
# ###############
# ##
# 
# # First re-label the factors
# TFPR.DATA$MA <- factor(TFPR.DATA$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), 
#                       labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
# TFPR.DATA$pooling <- factor(TFPR.DATA$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), 
#                       labels = ArtPOOLlabels)
# 
# # What is the point in number.thresholds most close to alpha = 0.05
# thresholds <- seq(0,1,length.out = number.thresholds)
# IDthreshold <- which(abs(c(thresholds - 0.05)) == min(abs(c(thresholds - 0.05))))
# 
# # Calculate the start and end of the line drops: most closest to IDthreshold
# lineStarts <- lineEnds <- c()
# for(i in 1:length(names_one)){
#   line <- TFPR.DATA[TFPR.DATA$source == names_one[i],'FP'][IDthreshold]
#   lineStarts <- c(lineStarts, line)
#   line_end <- TFPR.DATA[TFPR.DATA$source == names_one[i],'TP'][IDthreshold]
#   lineEnds <- c(lineEnds, line_end)
# }
# 
# # Add to data frame
# TFPR.DATA$LineStart <- rep(lineStarts, each = number.thresholds)
# TFPR.DATA$LineEnd <- rep(lineEnds, each = number.thresholds)
# 
# # Create all combinations of pooling and MA for labels of AUC
# combinations <- data.frame(expand.grid(
#   factor(poolmeth,levels = c('OLS', 'FixedEffect', 'MixedEffect'), 
#          labels = c('OLS', 'Fixed Effects', 'Mixed Effects')),
#   factor(metaMethods, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), 
#          labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))))
# colnames(combinations) <- c('pooling', 'MA')
# AUC_DAT <- data.frame(combinations,
#               AUC = format(round(unlist(AUC_VALUES),4), nsmall = 4))
# 



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










################################################################################
################################################################################
################################### APPENDIX ###################################
################################################################################
################################################################################

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


# Number of partial ROC figures in paper
partROCFigNumber <- c('A1','A2','A3')

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




