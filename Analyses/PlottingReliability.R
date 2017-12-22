####################
#### TITLE:     Calculate activation reliability of the MAs: plot in overlap matrix.
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

# Check if we are in the correct working directory
if(!grepl(pattern = 'Analyses', x = getwd())){
  stop('Working directory not at IntermediateResults')
}

# Location to save plots for the paper
LocFileSave <- '/Figures/'

# Specific functions from cowplot package.
# I gathered these from the cowplot Github page (https://github.com/wilkelab/cowplot).
# This is because the full package is not bug free on my machine as of 03/01/2017.
source('cowplot_functions.R')


##
###############
### Preparation
###############
##


# Libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(png)
library(lme4)
library(car)

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

# Dimension of the data
DIM <- c(53,63,46)

# Vector of pooling methods
poolmeth <- c(
	'FixedEffect',
	'OLS',
	'MixedEffect'
	)
numpoolmeths <- length(poolmeth)

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


##
###############
### Load values for overlap
###############
##

# data frame
# -- Overlap_K contains the values for the activation overlap
Overlap_K <- data.frame() %>% as_tibble()

for(f in 1:dim(DesignInfo)[1]){
  # WD gets f in loop 
  WD <- f
  
  # Select NRUNS, NSTUD and number of pairs
  NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
    unlist() %>% as.numeric()
  NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
    unlist() %>% as.numeric()
  
  # The pairwise comparissons
  combRuns <- c(1:NRUNS)
  PAIRS <- t(combn(combRuns,2))
  NPAIRS <- dim(PAIRS)[1]
  
  # Loop over the number of data frames to load in
  Overlap <- readRDS(paste('IntermediateResults/Overlap/Overlap_K_', NSTUD,'.rds', sep = ''))
  
  OverlapD <- data.frame('Overlap' = unlist(Overlap),
                         'pooling' = rep(poolmeth, NPAIRS * numMetaMethods),
                  'MA' = rep(names(Overlap), each = NPAIRS * numpoolmeths)) %>%
    `rownames<-`(1:c(NPAIRS * numMetaMethods * numpoolmeths)) %>%
    mutate(K = NSTUD, Pair = rep(1:NPAIRS, numMetaMethods * numpoolmeths))
  
  Overlap_K <- bind_rows(Overlap_K, OverlapD)
}

# Add labeling for plotting to factors
Overlap_K$pooling <- factor(Overlap_K$pooling, levels = sort(poolmeth), labels = sort(ArtPOOLlabels))
Overlap_K$MA <- factor(Overlap_K$MA, levels = sort(metaMethods), labels = sort(ArtLABMA))


##
###############
### Average overlap: plots and analyses
###############
##

Overlap_K %>%  
  ggplot(aes(x = K, y = Overlap, group = pooling)) +
  #geom_point(aes(colour = pooling)) + 
  geom_smooth(aes(colour = pooling), se = TRUE, method = 'loess') +
  facet_grid(~ MA) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  scale_color_brewer('Group level model', type = 'qual', palette = 2) +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))

Overlap_K %>% 
  ggplot(aes(x = factor(K), y = Overlap)) +
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

Overlap_K %>% group_by(pooling, MA, K) %>% 
  summarise(AvgOv = mean(Overlap)) %>%
  ggplot(aes(x = factor(K), y = AvgOv)) +
  geom_col(fill = 'gray', colour = 'black') +
  facet_grid(pooling ~ MA) +
  theme_bw(base_size = 9, base_family = "Helvetica") +
  theme(panel.grid.major = element_line(size = 0.2, color = "grey"),
        axis.line = element_line(size = 0.4, color = "black"),
        axis.title.y = element_text(size=11),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(size=10, angle = 45, vjust = 0.25, hjust = 0.5),
        strip.background = element_rect(colour = "white", fill = "white"))


# Fit mixed model with random intercept for pooling as this is clustered within MA
fitMM <- Overlap_K %>% lmer(Overlap ~ pooling * MA + K + (1|pooling), data = ., REML = FALSE)

# Anova with type II test
Anova(fitMM, type = 'II', test.statistic = 'Chisq')

# Adding other interaction terms
fitMMint <- Overlap_K %>% lmer(Overlap ~ pooling * MA * K + (1|pooling), data = ., REML = FALSE)
Anova(fitMMint, type = 'II', test.statistic = 'Chisq')

##
###############
### Figures in main text
###############
##

#### Overlap for K = 10, 20 and 35 ----
# Choose your WD to plot (in paper: 1, 6, 8)
WD <- 8

# Corresponding NRUNS and NSTUD
NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
  unlist() %>% as.numeric()
NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
  unlist() %>% as.numeric()

# Figure label in main text
OverlapFigNumber <- c(7,'','','','',8,'',9)

# The number of pairwise comparissons
combRuns <- c(1:NRUNS)
PAIRS <- t(combn(combRuns,2))
NPAIRS <- dim(PAIRS)[1]

# Load values of overlap
OverlapP <- readRDS(paste('IntermediateResults/Overlap/Overlap_K_', NSTUD,'.rds', sep = ''))

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
col1 <- c('#7fc97f','#beaed4','#fdc086')
col1 <- c('#08589e', '#7bccc4', '#ccebc5')

# Overlap heatmap
quartz(height = 6.5, width = 6.5,
	type = 'png', file = paste(getwd(), '/', LocFileSave, '/figure',OverlapFigNumber[WD],'_overlap_K',NSTUD,'.png', sep = ''), bg = 'white', canvas = 'white', dpi = 600)
ggplot(data = CORRM, aes(Var2, Var1, fill = value))+
geom_tile(color = "white") +
facet_grid(pooling ~ MA) +
geom_text(aes(label = str_replace(as.character(round(value,2)), "^0\\.", ".")), colour = "white", size = ifelse(WD == 1,2.5,4)) +
scale_fill_gradient2(low = col1[1], mid = col1[2], high = col1[3],
  midpoint = 0.5, limit = c(0,1), space = "Lab",
  name="Overlap\nCoefficient") +
	scale_x_discrete(name=paste("FOLD (when K = ", NSTUD, ")", sep=""),
			labels = 1:NRUNS) +
     scale_y_discrete(name="FOLD",
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



#### Average overlap ----

# Add predicted values to data frame: using mixed model fit
Overlap_K$Predict <- predict(fitMMint)
AvgOvData <- Overlap_K %>% group_by(pooling, MA, K, Predict) %>% 
  summarise(AvgOv = mean(Overlap)) %>% 
  gather(key = 'Type', value = 'Response', 4:5) %>% group_by(pooling, Type) %>%
  filter(Type == 'AvgOv')

# Plot average overlap with regression lines
quartz(height = 5, width = 8.2,
       type = 'png',file = paste(getwd(),'/', LocFileSave, 'figure10', 
                                 '_AVG_Overlap.png', sep = ''), dpi = 600)
Overlap_K %>% group_by(pooling, MA, K, Predict) %>% 
  summarise(AvgOv = mean(Overlap)) %>% 
  gather(key = 'Type', value = 'Response', 4:5) %>% group_by(pooling, Type) %>%
  ggplot(aes(x = K, y = Response)) +
  geom_line(aes(colour = pooling, linetype = Type, alpha = Type), size = 0.9) +
  geom_point(data = AvgOvData, aes(x = K, y = Response, 
                                   colour = pooling), size = 0.95) +
  facet_grid(~ MA) +
  scale_linetype_manual('Type', labels = c('Observed overlap', 'Fitted regression'),
                        values = c('solid', 'dashed')) +
  scale_x_continuous(name = 'Number of studies (K)') +
  scale_y_continuous(name = 'Average overlap of activation') +
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
### Figures in appendix
###############
##

#### Overlap for other K ----
# Choose your WD to plot (in paper: 1, 6, 8)
WDs <- c(2:5, 7)

for(a in 1:length(WDs)){
  WD <- WDs[a]
  # Corresponding NRUNS and NSTUD
  NRUNS <- DesignInfo %>% filter(WDs == WD) %>% select(FOLDS) %>% 
    unlist() %>% as.numeric()
  NSTUD <- DesignInfo %>% filter(WDs == WD) %>% select(K) %>% 
    unlist() %>% as.numeric()
  
  # Figure label in appendix
  OverlapFigNumber <- c('','A5','A5','A5','A5','','A5','')
  
  # The number of pairwise comparissons
  combRuns <- c(1:NRUNS)
  PAIRS <- t(combn(combRuns,2))
  NPAIRS <- dim(PAIRS)[1]
  
  # Load values of overlap
  OverlapP <- readRDS(paste('IntermediateResults/Overlap/Overlap_K_', NSTUD,'.rds', sep = ''))
  
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
  col1 <- c('#7fc97f','#beaed4','#fdc086')
  col1 <- c('#08589e', '#7bccc4', '#ccebc5')
  
  # Overlap heatmap
  quartz(height = 6.5, width = 6.5,
         type = 'png', file = paste(getwd(), '/', LocFileSave, 'Appendix/sup_',OverlapFigNumber[WD],'_overlap_K',NSTUD,'.png', sep = ''), bg = 'white', canvas = 'white', dpi = 600)
  print(ggplot(data = CORRM, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white") +
    facet_grid(pooling ~ MA) +
    geom_text(aes(label = str_replace(as.character(round(value,2)), "^0\\.", ".")), colour = "white", size = ifelse(WD == 1,2.5,4)) +
    scale_fill_gradient2(low = col1[1], mid = col1[2], high = col1[3],
                         midpoint = 0.5, limit = c(0,1), space = "Lab",
                         name="Overlap\nCoefficient") +
    scale_x_discrete(name=paste("FOLD (when K = ", NSTUD, ")", sep=""),
                     labels = 1:NRUNS) +
    scale_y_discrete(name="FOLD",
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
    coord_fixed())
  dev.off()
}

##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
