####################
#### TITLE:     File with all functions, to be sourced to
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/
#### First Modified: 01/07/2014
#### Notes:
#################


# NOTES:
# Contains list of functions to:
# - write local maxima in correct format for coordinate based meta-analysis
      # - Effect size based version (such as seed based d mapping)
      # - Activation Likelihood Estimation version
# - perform an effect size based meta-analyses (fixed and random effects models)
  #  - Estimation and inference!
# - function to apply transformations between voxel space and native space
# - function to perform voxelwise False Discovery Rate


# Some libraries
require(dplyr)
require(foreach)
require(doParallel)
require(compiler)
require(tcltk)
require(oro.nifti)
require(devtools)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STUDY LEVEL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


##
###############
### Pull local maxima from studies and write it to format for meta-analysis
###############
##

## Format for own Meta-Analysis
writeLocalMax <- function(numofsubj,indexofstudy,pathtofile,filename,pathtowrite,outputname,FirstColumn){
  ###------ HELP:
# numofsubj are the number of subjects in that particular study
# indexofstudy is the index (e.g. for the 8th study: index = 8)
# pathtofile is the absolute path to the file
# filename is how the input file is called
# pathtowrite is location where the file needs to be written
# outputname: how is output file called: WITHOUT THE .txt
# FirstColumn is the column with the x coordinate, y and z should be the next two columns

  if(length(grep(".txt",outputname))==1){stop("name of output file should not contain .txt!")}
  # Number of subjects
  S <- numofsubj
  # Index of study
  i <- indexofstudy
  # Open file
  dat <- read.table(paste(pathtofile,filename,sep=''),sep = "\t", header=TRUE, fill=TRUE,row.names=NULL)

  # Convert the Z-values to t-values: equating probabilities under the tails of the distributions (t'-> p->z')
    # Furthermore, we use the left side of the distribution to convert values at the extreme right side of the distribution,
    # since large Z-values result in Inf values (P-value becomes 1). However, a P-value of 0 is not possible. Hence we use the negative Z-value, calculate t-value and then
    # switch again to the correct sign.
  TVal <- (sign(dat[,2])*(-1))*qt(pnorm(-abs(dat[,2])),df=c(S-1))


  # Select the right columns and append amount of subjects next to it:
      # if there are no significant local maxima, then only write amount of subjects in study in the form of S S S S S S
  if(dim(dat)[1]==0){
    peak <- array(S,dim=c(1,6))
    colnames(peak) <- c("N","N","N","N","N","N")
    # Write to table
    write.table(peak, paste(pathtowrite,outputname,i,".txt",sep=""), sep = "\t", row.names=FALSE, col.names=TRUE, append=FALSE, quote=FALSE)
  } else {
    peak <- cbind(dat[,1],dat[,FirstColumn:(FirstColumn+2)],dat[,2],TVal,c(rep(S,dim(dat)[1])))
    names(peak) <- c("PeakIndex","X","Y","Z","ZVal","TVal","N")
    # Write to table
    write.table(peak[,c(2:7)], paste(pathtowrite,outputname,i,".txt",sep=""), sep = "\t", row.names=FALSE, col.names=TRUE, append=FALSE,quote=FALSE)
  }


}

# Format for ALE Meta-analysis
writeALELMax <- function(numofsubj,indexofstudy,pathtofile,filename,pathtowrite,outputname,individFile,FirstColumn){
  ###------ HELP:
# numofsubj are the number of subjects in that particular study
# indexofstudy is the index (e.g. for the 8th study: index = 8)
# pathtofile is the absolute path to the file
# filename is how the input file is called
# pathtowrite is location where the file needs to be written
# outputname: how is output file called: WITHOUT THE .txt
# individfile: logical if procedure needs to write the individal files of each study yes/no (TRUE/FALSE)
    # Mostly not really needed (ALE uses one file with blank line between studies).
    # WARNING: IF TRUE, THEN WE NEED A FOR LOOP FOR EACH STUDY (S)!
# FirstColumn specifies which is the first x column of the file (the next two should be y and z)

if(length(grep(".txt",outputname))==1){stop("name of output file should not contain .txt!")}
  # Number of subjects
  S <- numofsubj
  # Index of study
  i <- indexofstudy
  # Open file
  dat <- try(read.table(paste(pathtofile,filename,sep=''),sep = "\t", header=TRUE, fill=TRUE,row.names=NULL),silent=TRUE)

  # Write to table: in ALE format
        #// Reference=MNI
        #// HCP, 2014
        #X Y Z

  # Individual files (for each study) needed?
  if(individFile==TRUE){
    # First create the file
    peakFile <- file(paste(pathtowrite,outputname,i,".txt",sep=""))
    writeLines(c("// Reference=MNI","// HCP,2014",paste("// Subjects=",S,sep="")), peakFile)
    close(peakFile)
    # Then append the data
    if(dim(dat)[1]==0){
      # Include empty row
      cat("\n", file = paste(pathtowrite,outputname,i,".txt",sep=""))
    } else {
      write.table(dat[,FirstColumn:(FirstColumn+2)], paste(pathtowrite,outputname,i,".txt",sep=""), sep = "\t", row.names=FALSE, col.names=FALSE, append=TRUE )
    }
  }

  # Combine them in one file with all the coordinates of the studies
  # First iteration: create the file
  if(i==1){
    allPeaks <- file(paste(pathtowrite, "ALEPeaks.txt",sep=""))
    writeLines(c("// Reference=MNI","// HCP,2014_1",paste("// Subjects=",S,sep="")), allPeaks)
    close(allPeaks)
  }else{
    # Header for study
    cat("// Reference=MNI","\n",paste("// HCP,2014_",i,sep=""),"\n",paste("// Subjects=",S,sep=""),"\n",file = paste(pathtowrite, "ALEPeaks.txt",sep=""), sep="",append = TRUE)
    }
  # Append data in this master file: if there are local maxima!
    if(dim(dat)[1]==0){
      # Include empty row
      cat("\n", file = paste(pathtowrite, "ALEPeaks.txt",sep=""), append = TRUE)
    } else {
      write.table(dat[,FirstColumn:(FirstColumn+2)], paste(pathtowrite, "ALEPeaks.txt",sep=""), sep = "\t", row.names=FALSE, col.names=FALSE, append=TRUE )
      cat("\n", file = paste(pathtowrite, "ALEPeaks.txt",sep=""), append = TRUE)
    }

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PREPARING EFFECT SIZE BASED COORDINATE-BASED META-ANALYSIS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##
###############
### Helper functions for preparing the CBMA
###############
##


# Euclidean Distance: form = [x1,y1,z1,x2,y2,z2] to compare XYZ of voxel 1 with XYZ of voxel 2
ED <- function (LOC1,LOC2){
	x1 <- LOC1[1]
  	y1 <- LOC1[2]
  	z1 <- LOC1[3]
  	x2 <- LOC2[1]
  	y2 <- LOC2[2]
  	z2 <- LOC2[3]

    return((sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)))
}

# Unnormalized Gaussian Kernel K (Radua et al. 2011)
  KER <- function(D,FWHM){
  	return(exp(-(4*log(2))/(FWHM^2)*(D^2)))
  }


# Transform to MM distance:
	# transform xyz voxel distance to MM corresponding of the mm in the template (eg. MNI 152 2mm)
  	# e.g: MMDis(1,2) =
  	# euclidean distance of 1 voxel comes 2 mm in MNI 152 2mm template
  MMDis <- function(vox,dim){
  	return(vox*dim)
  }


# Masking
masking <- function(x,MASK){
  # X is data with NA values
  # Mask is binary with 0-1
  x[which(MASK==0)] <- NA
  return(x)
}


# Average weighting
  # Input should be: two vectors, one with the weights, the other with the to be averaged values (x)
  AvWe <- function(w,x){
    w[which(w==0)] <- 1
  	w[which(is.na(w))] <- 0
  	w <- w^2
	x[which(is.na(x))] <- 0

  	WProd <- w*x
  	WSum <- apply(WProd,c(1),sum)

  	WDen <- apply(w,c(1),sum)

  	weighted <- WSum/WDen
  	weighted[which(is.nan(weighted))] <- 0

	return(weighted)
  }

# XYZ <-- Voxel space
v <- 1
MNI <- c(45,63,36)
TAL <- c(41,59,32)
FSLVoxel <- c(0,0,0)
ifelse(v<3,origin<-MNI,origin<-TAL)

tovox <- function(x,origin){
  # Choose between MNI, TAL and FSLVox.
    # FSLVox is just voxelspace in FSL but it starts at 0,0,0 while in R at 1,1,1, hence everything needs to be shifted one voxel (x+1,y+1,z+1)
    if(sum(origin)==0){
      c(x+1)
    }else{
    (trunc(c(-1,1,1)*x/2+origin))
    }
  }

# Voxel space <-- XYZ space
toxyz<-function(x,origin) ((x-origin)*2)*c(-1,1,1)

# Load in new studies
loadNewStud <- function(S,loc,name){
  allStud <- list()
  for(i in 1:S){
    stud <- read.csv(paste(loc,name,i,'.txt',sep=''),sep="",header=TRUE)
    allStud <- c(allStud, i = list(stud))
  }
  return(allStud)
}



####################################
# Global function to prepare studies
####################################


 prepare <- function(STUD, DIM, templateMM, voxelSize, FWHM, MASK, origin=MNI, transformation="tovox"){
    Num.S <- dim(summary(STUD))[1]

    brain <- array(NA,dim=c(prod(DIM)))

    DiMap <- expand.grid(1:DIM[1],1:DIM[2],1:DIM[3])


    # For loop: for S studies
    for(i in 1:Num.S){
      print(paste("At Study ", i, " of ", Num.S,sep=""))
      Foc.Study <- STUD[[i]]
        Foci.Study <- Foc.Study
          if(colnames(Foci.Study)[1]=="N"){
            WeES.Study <- array(0,dim=prod(DIM))
          } else {
            if(transformation=="tovox"){
                Foci.Study[,c('X','Y','Z')] <- t(apply(Foc.Study[,c('X','Y','Z')],c(1),tovox,origin=origin))    # Take the foci of each study and convert them to voxel space
              } else if(transformation=="IMAGEN") {
                tmp.Foci.Study <- apply(t(Foc.Study[,c('X','Y','Z')]),c(2),MNItoR,TM=TransMatrix)               # Take the foci of each study and convert them to voxel space
                Foci.Study[,c('X','Y','Z')] <- t(tmp.Foci.Study[-4,])
              }
      dim.Study <- dim(Foci.Study)[1]

      DiMap.Study <- list()
        for(j in 1:dim.Study){
          DiMap.Study <- c(DiMap.Study, j = list(DiMap))
          }

      LisFocStu <- as.list(data.frame(t(Foci.Study[,c('X','Y','Z')])))
      tvals <- as.list(Foci.Study$TVal)

      EUDis <- mapply(ED, DiMap.Study, LisFocStu)

      TemMM <- MMDis(EUDis, templateMM)
        TemMM[which(TemMM>20)] <- NA

      KERNEL <- lapply(data.frame(TemMM), KER, FWHM=FWHM)

      MASK <- array(MASK,dim=prod(dim(MASK)))
      KERNEL.MASKED <- lapply(KERNEL,masking,MASK=MASK)

      HeG <- as.list(mapply(hedgeG, tvals,Foc.Study[1,'N']))

      # Here we need to switch to for loop as studies with too much foci consume too much memory for listwise calculations...
      ES.Study <- c()
        for(s in 1:length(KERNEL.MASKED)){
          MASKED.TMP <- KERNEL.MASKED[[s]]
          DIST <- MASKED.TMP * HeG[[s]]
          ES.Study <- c(ES.Study, DIST)
        }
      ES.Study <- array(ES.Study,dim=c(prod(DIM),length(KERNEL.MASKED)))

      WeES.Study <- AvWe(TemMM,ES.Study)
      }

      brain <- cbind(brain,WeES.Study)
    }
    brain <- brain[,-1]
    return(brain)
  }




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# EXECUTING META-ANALYSIS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


##
###############
### Helper functions for ES CBMA meta-analysis
###############
##

# Calculate Hedges'g from t-statistic from one sample t-test with N subjects
hedgeG <- function(t,N){
  J <- 1-(3/((4*(N-1))-1))
  G <- (t/sqrt(N))*J
  return(G)
}

varHedge <- function(g,N){
  value <- (1/N) + (1 - (gamma((N - 2) / 2) / gamma((N - 1) / 2))^2 * (N - 3) / 2) * g^2
    return(round(value,7))
}

# Calculate between-studies variance
tau <- function(Y,W,k){
  C <- sum(W)-(sum(W^2)/sum(W))
  df <- k-1
  Q <- sum(W*Y^2)-sum(W*Y)^2/sum(W)
  if(Q < df){
    T2 <- 0
  }else{
    T2 <- (Q-df)/C
  }
  return (T2)
}

# Calculate weighted mean
wMean <- function(Y,W){
	M <- sum(W*Y)/sum(W)
	return(M)
}


##################################################################
# Global function to perform meta-analysis: Mixed Effects Analysis
##################################################################

metaAnMix <- function(allStud, brain, DIM){
  # Extract amount of participants in each study
  N.S <- lapply(allStud, "[", 1, 'N')

  # Calculate variance and then weights
  VarWS <- mapply(varHedge, as.list(as.data.frame(brain)),N.S)
  weigFix <- 1/VarWS

  # Calculate the between subject variance: tau.
    HeGL <- as.list(as.data.frame(t(brain)))
    weigFixL <- as.list(as.data.frame(t(weigFix)))
    K <- length(N.S)
  VarBS <- as.vector(mapply(tau,Y=HeGL,W=weigFixL,k=K))

  # Random effects weights:1/(VarWS+VarBS)
  weigRan <- 1/(apply(VarWS,c(2),function(x,y){x+y},y=VarBS))
    weigRanL <- as.list(as.data.frame(t(weigRan)))
  # Weighted mean
  WMeanRan <- as.vector(mapply(wMean,HeGL,weigRanL))

  # Return to dimension of template
  WMeanRan[which(is.na(WMeanRan))] <- 0
  WMeanRan <- array(WMeanRan,dim=DIM)

  # Variance of summary effect: reciprocal of sum of weights
  varRan <- 1/apply(weigRan, c(1), sum)
  seRan <- sqrt(varRan)

  # Z-statistic
  ZstatRan <- WMeanRan / seRan

  combined <- list(WMeanRan,varRan,VarBS, ZstatRan)

  return(combined)
}



############################################################################
### Generic function to perform meta-analysis through FIXED Effects Analysis
############################################################################

metaAnFix <- function(allStud, brain, DIM){
  # Extract amount of participants in each study
  N.S <- lapply(allStud, "[", 1, 'N')

  # Calculate variance and then weights
  VarWS <- mapply(varHedge, as.list(as.data.frame(brain)),N.S)
  weigFix <- 1/VarWS

    # Puth the ES and weights in a list
    HeGL <- as.list(as.data.frame(t(brain)))
    weigFixL <- as.list(as.data.frame(t(weigFix)))

  # Weighted mean
  WMeanFix <- as.vector(mapply(wMean,HeGL,weigFixL))

  # Return to dimension of template
  WMeanFix[which(is.na(WMeanFix))] <- 0
  WMeanFix <- array(WMeanFix,dim=DIM)

  # Variance of summary effect: reciprocal of sum of weights
  varFix <- 1/apply(weigFix, c(1), sum)
  seFix <- sqrt(varFix)

  # Z-statistic
  ZstatFix <- WMeanFix / seFix

  combined <- list(WMeanFix, ZstatFix)

  return(combined)
}



#####################################################################
# Functions to calculate null distribution permutation based P-values
#####################################################################


# Function to calculate null distribution based on permutations: TYPE: 1 == FIXED, 2 == MIXED Effects
MakePermDis <- function(brain,allStud,type,mask){
  TRANSmask <- array(mask,dim=c(prod(dim(mask)),1))
  maskedIND <- which(TRANSmask==1,arr.ind=TRUE)[,1]
  ITER <- length(maskedIND)
  PermBrain <- array(NA,dim=c(length(maskedIND),dim(brain)[2]))
    for(i in 1:ITER){
      indices <- sample(maskedIND,dim(brain)[2], replace=TRUE)
      PermBrain[i,] <- diag(brain[indices,c(1:dim(brain)[2])])
    }
  if(type==1){
      PermDat.tmp <- metaAnFix(allStud,PermBrain,DIM=ITER)[[2]]
      PermDat <- array(PermDat.tmp,dim=ITER)
    }else{
      PermDat.tmp <- metaAnMix(allStud,PermBrain,DIM=ITER)[[4]]
      PermDat <- array(PermDat.tmp,dim=ITER)
    }

    return(PermDat)
  }
# Make it precompiled
CMakePermDis <- cmpfun(MakePermDis)



##################################################################################
# Functions to calculate permutation based P-values, need permutation distribution
##################################################################################

# Remove abundant zeroes and replace it with just one (for permutations)
One0 <- function(x){c(x[which(x>0)],0)}

# Function to have the permutations in shorter vector this is: using only the unique values
  # so it does not have to do the computations over and over
PreparePerm <- function(permDistr){
  PermOne0 <- One0(permDistr)       # Remove all abundant zeroes
  Ra <- rank(PermOne0,ties="min")     # Define the ranks of the sorted permutations
  UnRa.tmp <- unique(Ra)          # use only the unique values
  UnRa <- max(Ra)-UnRa.tmp + 1      # Reverse rank (so highest mean -> lowest p-value)
  UnPerm <- unique(PermOne0)        # Use unique permutation based values

  TotLen <- length(permDistr)       # Total length (because we remove so much zeroes)
  return(list(UnPerm,UnRa,TotLen))
}

# Actual function to calculate permutation based p-value
  # Need to decide whether an observed value > max(permutations) should have p-value = 0 or 1/length(Perm)
PermPValRank <- function(Emp,UnPerm,UnRa,TotLen){
  matchEmPer <- which(abs(Emp-UnPerm) == min(abs(Emp-UnPerm)),arr.ind=TRUE)
  matchRank <- UnRa[matchEmPer][1]

  Pval <- ifelse(matchRank==max(UnRa),1,matchRank/TotLen)
  return(Pval)
}

# Same function, but pre-compiled. To speed it up.
CPermPValRank <- cmpfun(PermPValRank)

# Improved version of the function:
  # Using counts of the null distribution
  # Better accuracy, same speed advantage!
  # Input:
      # - Emperival values without zero (these get P = 1 anyway)
      # - FULL permutation based null-distribution
PermPValCountVectorized <- function(EmpWMeanExNull, PermWMean){
  # First we create a function that sums the counts in a given null distribution, larger than or equal than the observed, emperical value.
  CountFU <- function(EmpValue, CountedPerm){
    # Take the ID of the values larger than or equal than emperical value
    LargIDs <- EmpValue <= CountedPerm[,1]
    # Sum the counts of the selected elements in the distribution
    PValCountFU <- (sum(CountedPerm[LargIDs,2]) + 1) / (length(PermWMean) + 1)
    return(PValCountFU)
  }
  # Now for real, count the provided null distribution.
  CountedPerm <- count(PermWMean)
  # Apply the CountFU function from above to the emperical values, using the counted permutation distribution.
  PValExNullCountVect <- apply(EmpWMeanExNull, 1, CountFU, CountedPerm = CountedPerm)
  return(PValExNullCountVect)
}

# Same function, but pre-compiled. To speed it up.
CPermPValCountVectorized <- cmpfun(PermPValCountVectorized)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# VARIA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##
###############
### Remove values from a pool of sampled values
###############
##

remove <- function(keep,remove){
  for(i in 1:length(remove)){
    keep <- keep[which(keep != remove[i])]
  }
  return(keep)
}


##
###############
### Converting IMAGEN data from FSL/or when using R to MNI format
###############
##

#  We have two transformation functions :
#           the first one is for using coordinates from an SPM (coordinates e.g. from FSL analysis go from 0 -> i-1).
#               ==> use SPMtoMNI function
#           the second RtoMNI is for using R coordinates in an analysis (because R indexing goes from 1 -> i)
#               ==> use RtoMNI function
#  Both should use the TransMatrix, which is the qto_xyz transformation matrix
TransMatrix <- matrix(c(-3,0,0,78,0,3,0,-112,0,0,3,-50,0,0,0,1),byrow=TRUE,ncol=4)

# SPM to MNI:
SPMtoMNI <- function(TM,co){
  matMult <- TM%*%c(co,1)
  cor <- matMult+(matMult%%3)/3
  return(trunc(cor))
  }
# R to MNI: subtracts one value from the coordinates because R goes from 1 -> i while FSL goes from 0 -> i-1
RtoMNI <- function(TM,co){
  matMult <- TM%*%c(co-1,1)
  cor <- matMult+(matMult%%3)/3
  return(trunc(cor))
  }

# Approximation to exact transformation
floorSPMtoMNI <- function(TM,co){floor(TM%*%c(co,1)+0.5)}

  # Usage of these matrix with function:
  # Say we have coordinates (10 15 20) and (30 15 0), from FSL cluster analysis.
    # Put these in matrix with in the columns (!) the coordinates:
    # coord <- matrix(c(10,15,20,30,15,0), ncol = 2)
    # Then use apply over the columns: TM is the transformation matrix and co are the coordinates
    # apply(coord, 2, SPMtoMNI, TM=TransMatrix)

# Functions with inverse of transformation matrix, if one wants to go back
MNItoSPM <- function(TM,co){
  TM <- solve(TM)
  matMult <- TM%*%c(co,1)
  cor <- matMult+(matMult%%3)/3
  return(trunc(cor))
  }

# For usage in R
MNItoR <- function(TM,co){
  TM <- solve(TM)
  matMult <- TM%*%c(co,1)
  cor <- matMult+(matMult%%3)/3
  return(trunc(cor+1))
  }


##
###############
### Calculating the FDR with either assumption of independence/positive dependence or no assumption
###############
##

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


# Copy header from source to target nifti image
CopyHeader <- function(target, source, pixdim = c(-1,3,3,3,1,1,1,1)){
	EndNIFTI <- nifti(img = target,
				sizeof_hdr = source@sizeof_hdr,
				datatype = source@datatype,
				db_name = source@db_name,
				extents = source@extents,
				regular = source@regular,
				dim_info = source@dim_info,
				dim_ = source@dim_,
				intent_p1 = source@intent_p1,
				intent_p2 = source@intent_p2,
				intent_p3 = source@intent_p3,
				intent_code = source@intent_code,
				bitpix = source@bitpix,
				slice_start = source@slice_start,
				pixdim = source@pixdim,
				vox_offset = source@vox_offset,
				scl_slope = source@scl_slope,
				scl_inter = source@scl_inter,
				slice_end = source@slice_end,
				slice_code = source@slice_code,
				xyzt_units = source@xyzt_units,
				slice_duration = source@slice_duration,
				toffset = source@toffset,
				glmax = source@glmax,
				glmin = source@glmin,
				descrip = source@descrip,
				aux_file = source@aux_file,
				sform_code = source@sform_code,
				quatern_b = source@quatern_b,
				quatern_c = source@quatern_c,
				quatern_d = source@quatern_d,
				qoffset_x = source@qoffset_x,
				qoffset_y = source@qoffset_y,
				qoffset_z = source@qoffset_z,
				srow_x = source@srow_x,
				srow_y = source@srow_y,
				srow_z = source@srow_z,
				intent_name = source@intent_name,
				magic = source@magic,
				extender = source@extender,
				reoriented = source@reoriented
				)
			# For some reason, the qform_code cannot be written directly. We add it here
			EndNIFTI@qform_code <- source@qform_code
			EndNIFTI@pixdim <- pixdim
	return(EndNIFTI)
}




