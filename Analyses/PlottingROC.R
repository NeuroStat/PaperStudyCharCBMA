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

# Location to save plots for the paper
LocFileSave <- '/Figures/'

# Libraries
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(lme4)
library(car)

##
###############
### Preparation
###############
##


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

# 3 data frames:
# -- TFPR.DATA_K contains the true and false positive rates
# -- AUC_avg is the Area Under the Curve of the average ROC curves
# -- AUC_long is the Area Under the Curve of each ROC curve
TFPR.DATA_K <- AUC_avg <- AUC_long <- data.frame() %>% as_tibble()

for(f in 1:dim(DesignInfo)[1]){
  # WD gets f in loop 
  WD <- f
  
  # Select NRUNS, NSTUD
  NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
    unlist() %>% as.numeric()
  NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
    unlist() %>% as.numeric()
  
  # Loop over the number of data frames to load in
  TFPR.DATA <- readRDS(paste('IntermediateResults/TFPR/TFPR_K_', NSTUD,'.rds', sep = ''))
  AUC_VALUES_avg <- readRDS(paste('IntermediateResults/AUC_avgFPTP/AUC_VALUES_K_', NSTUD,'.rds', sep = ''))
  AUC_VALUES_long <- readRDS(paste('IntermediateResults/AUC_long/AUC_long_K_', NSTUD,'.rds', sep = ''))
  
  # First re-label the factors
  TFPR.DATA$MA <- factor(TFPR.DATA$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), 
                         labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
  TFPR.DATA$pooling <- factor(TFPR.DATA$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), 
                              labels = ArtPOOLlabels)
  AUC_VALUES_long$MA <- factor(AUC_VALUES_long$MA, levels = c('ALE', 'FixedEffUn', 'RanEffUn'), 
                         labels = c('ALE', 'Fixed Effects MA','Random Effects MA'))
  AUC_VALUES_long$pooling <- factor(AUC_VALUES_long$pooling, levels = c('OLS', 'FixedEffect', 'MixedEffect'), 
                              labels = ArtPOOLlabels)
  
  # What is the point in number.thresholds most close to alpha = 0.05, needed for line drops
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
                AUC = AUC_VALUES_avg %>% unlist() %>% 
                    round(.,4) %>% format(., nsmall = 4) %>% as.numeric())
  
  # Add K to data frame and convert to tibble
  AUC_DAT <- AUC_DAT %>% mutate(K = NSTUD) %>% as_tibble()

  # Bind TFPR and both versions of AUC to data frames
  TFPR.DATA_K <- bind_rows(TFPR.DATA_K, TFPR.DATA)
  AUC_avg <- bind_rows(AUC_avg, AUC_DAT)
  AUC_long <- bind_rows(AUC_long, as_tibble(AUC_VALUES_long))
  
  # Remove objects
  rm(TFPR.DATA, AUC_VALUES_avg, AUC_DAT, AUC_VALUES_long)
}


##
###############
### Analysis on AUC: several plots
###############
##

AUC_long %>% select(-Source) %>% 
  ggplot(aes(x = K, y = AUC, group = pooling)) +
  #geom_point(aes(colour = pooling)) + 
  geom_smooth(aes(colour = pooling), se = TRUE, method = 'loess') +
  facet_grid(~ MA) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))

AUC_long %>% select(-Source) %>% 
  ggplot(aes(x = K, y = AUC, group = pooling)) +
  geom_point(aes(colour = pooling), size = 0.75, alpha = 0.75) + 
  geom_smooth(aes(colour = pooling), se = FALSE, method = 'loess') +
  facet_grid(~ MA) +
  scale_color_brewer('Group level model', type = 'qual', palette = 6) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))


AUC_long %>% select(-Source) %>% 
  ggplot(aes(x = factor(K), y = AUC)) +
  geom_boxplot() +
  facet_grid(pooling ~ MA) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))


AUC_long %>% select(-Source) %>% group_by(pooling, MA, K) %>% 
  summarise(AvgAUC = mean(AUC)) %>%
  ggplot(aes(x = factor(K), y = AvgAUC)) +
  geom_col() +
  facet_grid(pooling ~ MA) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))

AUC_long %>% group_by(pooling, MA) %>%
  ggplot(aes(x = K, y = AUC)) +
  geom_point(aes(colour = pooling, shape = MA), size = 0.9) +
  scale_x_continuous(name = 'Number of studies (K)') +
  scale_y_continuous(name = 'AUC') +
  scale_shape_manual('Model for CBMA', values = c(1, 15, 20)) +
  scale_color_brewer('Group level model', type = 'qual', palette = 6) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))


AUC_long %>% group_by(pooling, MA, K) %>% 
  summarise(AvgAUC = mean(AUC)) %>%
  ggplot(aes(x = K, y = AvgAUC)) +
  geom_line(aes(colour = pooling, linetype = MA), size = 0.9) +
  scale_x_continuous(name = 'Number of studies (K)') +
  scale_y_continuous(name = 'AUC') +
  scale_linetype_manual('Model for CBMA', values = c('solid', 'dotted', 'dashed')) +
  scale_color_brewer('Group level model', type = 'qual', palette = 6) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))


AUC_long %>% group_by(pooling, MA, K) %>% 
  summarise(AvgAUC = mean(AUC)) %>%
  ggplot(aes(x = K, y = AvgAUC, group = pooling)) +
  geom_line(aes(colour = pooling), size = 0.9) +
  facet_grid(~ MA) +
  scale_x_continuous(name = 'Number of studies (K)') +
  scale_y_continuous(name = 'Area Under the Curve (AUC)') +
  scale_color_brewer('Group level model', type = 'qual', palette = 6) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))


# Fit mixed model with random intercept for pooling as this is clustered within MA
fitMM <- lmer(AUC ~ pooling * MA + K + (1|pooling), data = AUC_long, REML = FALSE)
# One without random intercept
fitLM <- lm(AUC ~ pooling * MA + K, data = AUC_long)

# AIC is better for regular LM
AIC(fitLM)
AIC(fitMM)

# Anova with type II test
Anova(fitLM, type = 'II', test.statistic = 'Chisq')
Anova(fitMM, type = 'II', test.statistic = 'Chisq')

# Other interaction effects
fitMMint <- lmer(AUC ~ pooling * MA * K + (1|pooling), data = AUC_long, REML = FALSE)
Anova(fitMMint, type = 'II', test.statistic = 'Chisq')

# Some summary values
AvgAUCData %>% group_by(pooling, MA) %>% 
  summarise(OverK = mean(Response))
AvgAUCData %>% group_by(MA, K) %>% 
  summarise(OverGM = mean(Response))


##
###############
### Figures in main text
###############
##

#### ROC for K = 10, 20 and 35 ----
# Choose your WD to plot (in paper: 1, 6, 8)
WD <- 8

# Corresponding NSTUD
NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
  unlist() %>% as.numeric()

# Number of ROC figures in paper
ROCFigNumber <- c(3,'','','','',4,'',5)

# Objects to plot
TFPR.DATA <- TFPR.DATA_K %>% filter(K == NSTUD)
AUC_DAT <- AUC_avg %>% filter(K == NSTUD)

quartz(height = 6.5, width = 3.5,
	type = 'png',file = paste(getwd(),'/', LocFileSave, 'figure', 
	       ROCFigNumber[WD],'_ROC_K',NSTUD,'.png', sep = ''), dpi = 600)
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


#### Average ROC's ----

# Add predicted values to data frame: using mixed model fit
AUC_long$Predict <- predict(fitMM)
AvgAUCData <- AUC_long %>% group_by(pooling, MA, K, Predict) %>% 
  summarise(AvgAUC = mean(AUC)) %>% 
  gather(key = 'Type', value = 'Response', 4:5) %>% group_by(pooling, Type) %>%
  filter(Type == 'AvgAUC')

# With fitted regression lines
quartz(height = 5, width = 8.2,
       type = 'png',file = paste(getwd(),'/', LocFileSave, 'figure6', 
                                 '_AVG_ROC.png', sep = ''), dpi = 600)
AUC_long %>% group_by(pooling, MA, K, Predict) %>% 
  summarise(AvgAUC = mean(AUC)) %>% 
  gather(key = 'Type', value = 'Response', 4:5) %>% group_by(pooling, Type) %>%
  ggplot(aes(x = K, y = Response)) +
  geom_line(aes(colour = pooling, linetype = Type, alpha = Type), size = 0.9) +
  geom_point(data = AvgAUCData, aes(x = K, y = Response, 
                                    colour = pooling), size = 0.95) +
  facet_grid(~ MA) +
  scale_linetype_manual('Type', labels = c('Observed AUC', 'Fitted regression'),
                        values = c('solid', 'dashed')) +
  scale_x_continuous(name = 'Number of studies (K)') +
  scale_y_continuous(name = 'Area Under the Curve (AUC)', limits = c(0.79,0.907)) +
  scale_color_brewer('Group level model', type = 'qual', palette = 2) +
  scale_alpha_manual('', values = c(1,0.5), guide = 'none') + 
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))
dev.off()

##
###############
### Figure with all ROC's for appendix
###############
##

ggplot(TFPR.DATA_K, aes(FP,TP, group = K)) +
  geom_point(aes(colour = K), size = 0.05, alpha = 0.7) +
  facet_wrap(MA ~ pooling, dir = "v") +
  scale_x_continuous(name = paste('False Positive Rate', sep = '')) +
  scale_y_continuous(name = paste('True Positive Rate', sep = '')) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  ggtitle('Average ROC curves for each K') +
  theme(panel.grid.major = element_line(size = 0.3, color = "#bdbdbd"),
        panel.grid.minor = element_line(size = 0.3, color = "#f0f0f0"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=9.2),
        axis.title.x = element_text(size=9.2),
        axis.text.x = element_text(angle = 45, vjust = 0.25,hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))

ggsave(file = paste(getwd(), LocFileSave, 'Appendix/sup_A1_ROC_all_K.pdf', sep = ''), device = "pdf",
       width = 6, height = 7, units = 'in', dpi = 300)





