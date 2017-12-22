#!/bin/sh

####################
#### TITLE:     -------------Continuous K STUDIES VERSION------------
####			Program to sample subjects over all studies from the Imagen data pool, pool them through fixed effects, ols or mixed effects.
####			And then perform the fixed and random effects meta-analysis (ALE in seperate scripts).
####			This is repeated for several runs/folds.
#### Author: Han Bossier
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_FixRan/FixRanStudy/Imagen/35Studies
#### First Modified: 20/11/2017 (ContK version, based on original of 13/11/2014)
####
#### Notes: some terminology is outdated and got mixed (eg. runs = folds)
#################

#### Structure of folders:
# Scripts (location of the main script and all scripts it uses)
#--# Working Directory (called Results were studies (after pooling), meta analyses and group analysis are stored in each time a folder for a run/fold)
#---# The number of studies in the MA (K), written as K_...
#-----# The number of RUN/FOLD we are in (e.g. Run_1)
#--------# Studies
#------------# e.g. Study_1
#----------------# FixedEffect
#----------------# OLS
#----------------# MixedEffect
#---------# MetaAnalyses
#------------# e.g. FixedEff
#----------------# FixedEffect
#----------------# OLS
#----------------# MixedEffect
#---------# Group Analysis

# Actual data is in other location!


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pre Step: Global variables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	# Array of amount of subjects in each study
	arrSub=()
	# Which computer: HPC or MAC
	COMPUTER=HPC
	# Location of the data
	if [ "$COMPUTER" = MAC ] ; then
		data=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/Renamed/processed/data
	fi
	if [ "$COMPUTER" = HPC ] ; then
		# Which is your vsc number
		vsc=${13}
		data=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/Imagen/data
	fi
	# Location of mask for the meta-analysis
	locMask=($data/AllSubMask.nii)
	# How many studies are there in the MA? (K or NS, number of studies)
	NS=${14}
	# Number of subjects in the group analysis
	groupNS=${15}
	# Pooling methods considered and number of pooling methods
	PoolMeth=(FixedEffect OLS MixedEffect)
	NumPoolMet=${#PoolMeth[@]}
		# We need to subtract 1 because we will start our for loops at zero (arrays in bash are defined starting from 0, 1, 2,...)
		NumPoolMet=$(($NumPoolMet - 1))
	# Which thresholding method for group analysis?
	threshGroup=FDR
	# Significance level in group analysis
	signLevel=0.001
	# Which thresholding method for individual studies?
	threshStudy=FDR
	# Excursion set for peak level inference in individual studies
	excursionSet=2
	# Significance level for individual studies.
	signLevelStud=0.05
	# Minimal cluster size for individual studies.
	minclustersizeStud=1
	# Which template do you want to use in the ALE program?
	ALETEMPLATE='MNI_gm_10'
	# Does MATLAB version of ALE work on HPC?
	MATLABVERSION=NO
	# In which resampling run (FOLD) are we?
	RUN=${1}
	# Location of the identifiers for subject and scanner location
	if [ "$COMPUTER" = MAC ] ; then
		AllIDs=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/OurIDs
	fi
	if [ "$COMPUTER" = HPC ] ; then
		AllIDs=$data/OurIDs
	fi




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step one: Preparation
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Location of scripts folder
if [ "$COMPUTER" = MAC ] ; then
	SCRPT=/Users/hanbossier/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Imagen/ContK
fi
if [ "$COMPUTER" = HPC ] ; then
	SCRPT=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/ContKStudiesFast
fi
cd "${SCRPT}"

### Location of results
if [ "$COMPUTER" = MAC ] ; then
	# Create folder for run (formerly called Results, but now on external HD)
	RUNWD=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/10RunsSites/CognitiveTaskContK
fi
if [ "$COMPUTER" = HPC ] ; then
	# Check if Results already exists, if not, create
	if [ ! -d "Results" ]; then
	  mkdir Results
	fi
	# Same for K_$NS folder
	if [ ! -d "Results/K_"$NS ]; then
		mkdir Results/K_$NS
	fi
	RUNWD=$SCRPT/Results/K_$NS
fi


# Create folder for that specific run
cd "${RUNWD}"
mkdir Run_$RUN
# Define workspace
WD=$RUNWD/Run_$RUN
cd "${WD}"

# Make and define Studies, MetaAnalyses, GroupAnalysis and SubMetGroupAn folder (the latter one is to compare
#	group study of subjects used in the meta-analysis with the GT (which is other half of subjects))
mkdir Studies
mkdir MetaAnalyses
mkdir GroupAnalysis
mkdir SubMetGroupAn
Studies=$WD/Studies
MetaAnalyses=$WD/MetaAnalyses
GroupAnalysis=$WD/GroupAnalysis
SubMetGroupAn=$WD/SubMetGroupAn


# Need to have a temporary folder, one in each run, so there is no overlap in calculations (parallel processing)!
	# Create temp folder
	mkdir $WD/temp
	temp=$WD/temp




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step two: make the directories corresponding to each study, then go into the folder and make
# 			fixed effects, OLS, mixed effects pooling
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd "${Studies}"

# If S2FOLDERS is TRUE then execute (jumps to fi)
if [ ${2} = TRUE ] ; then

for i in $(eval echo "{1..$NS}"); do
 mkdir Study_$i
 cd Study_$i
 mkdir FixedEffect
 mkdir OLS
 mkdir MixedEffect
 cd ..
done

# S2FOLDERS
fi



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step three: Sampling
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S3SAMPLING is TRUE then execute (jumps to fi)
if [ ${3} = TRUE ] ; then

echo Sampling subjects
cd "${SCRPT}"
Rscript "${SCRPT}"/LoadSubjContK.R "$Studies" "$RUN" "${SCRPT}"/_K_"${NS}" "${NS}" &> ROutput_SamplingCrossVal.txt


# S3SAMPLING
fi




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step four: pooling subjects
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Always execute this part of the code

# Subglobal variables
subjects=0
flameoMethod=(fe ols flame12)

# For loop over all studies to make an array of amount of subjects in each study:
for i in $(eval echo "{1..$NS}"); do
	cd "${Studies}"
	echo Getting subjects array:....... $i
	cd Study_$i/

	IFS=$'\n' read -d '' -r -a subjects < study_$i.txt

	NumSub=${#subjects[@]}
	arrSub=("${arrSub[@]}" $NumSub)			# Array of subjects, used for the meta-analysis (do not unlist!)

	# Re-initialize variable subject (not the arrSub variable)
	unset -v 'subjects'
done

# Go back to studies directory
cd "${Studies}"

# Print amount of subjects
echo "Amount of subjects in all studies = '${arrSub[@]}'"



# If S4POOL is TRUE then execute (jumps to fi)
if [ ${4} = TRUE ] ; then

# For loop over all studies and then each method of pooling the subjects
for i in $(eval echo "{1..$NS}"); do
	echo Pooling study:....... $i
	cd Study_$i/
	WRITE=$Studies/Study_$i/
	ConFile=
	VarConFile=

	# Get subject ID's
	IFS=$'\n' read -d '' -r -a subjects < study_$i.txt
	# ID variable ONLY for looping purpose (if number of subjects = 12: then BASH needs to loop from 0-11 !!!)
	NumSub=${#subjects[@]}
	NumSub=$(($NumSub - 1))

	# Go to the temp folder. Copy the data from the subjects to it.
	cd "${temp}"
	for j in $(eval echo "{0..$NumSub}"); do
		SubID=${subjects[$j]}
		scp -r $data/$SubID/SessionB/EPI_global/swea/con_0023_$SubID.nii .
		scp -r $data/$SubID/SessionB/EPI_global/swea/varcon_0023_$SubID.nii .
		# Copy the universal mask to the temp folder, rename it with a subject index
		scp -r $data/AllSubMask.nii .
		mv AllSubMask.nii AllSubMask_$SubID.nii
		# Arrays
		ConFile=("${ConFile[@]}" con_0023_$SubID.nii)
		VarConFile=("${VarConFile[@]}" varcon_0023_$SubID.nii)
	done

	# For loop over all methods considered.
		# Use fslmerge with the corresponding method and then extract the files.
		for k in $(eval echo "{0..$NumPoolMet}"); do
			MET=${PoolMeth[$k]}
			OUTcope=$Studies/Study_$i/"$MET"/cope
			OUTvarcope=$Studies/Study_$i/"$MET"/varcope
			# Move mask to here
			scp -r "${temp}"/AllSubMask_* "${Studies}"/Study_$i/"$MET"/
			# Continue with merging the 3D images
			$FSLDIR/bin/fslmerge -t "$OUTcope" ${ConFile[@]}
			$FSLDIR/bin/fslmerge -t "$OUTvarcope" ${VarConFile[@]}
			cd "${WRITE}${MET}"
			gunzip cope.nii.gz
			gunzip varcope.nii.gz
			# Now create the needed files for flameo (design.mat, design.con and design.grp and a mask) through R file designCrossVal.R
			cd "${Studies}"
			Rscript "${SCRPT}"/designCrossVal.R "$WRITE$MET/" "$(($NumSub + 1))" &> ROutput_designCrossVal.txt

			# Merge the masks: 4D mask
			cd "${WRITE}${MET}"
			$FSLDIR/bin/fslmerge -t "${Studies}"/Study_$i/"$MET"/mask.nii AllSubMask_*.nii

			# Remove unzipped copes, if they are present (flameo doesn't want two sets of copes in the same folder)
			if [ -f "cope.nii.gz" ]
				then
        		rm -r cope.nii.gz
			fi

			# Remove individual masks
			rm -r AllSubMask_*.nii

			# Now perform flameo to pool subjects: cope is a 4D data structure (4th dimension are the subjects)
			echo .......pooling.........
			if [ "$MET" = OLS ] ; then
				$FSLDIR/bin/flameo --cope=cope --mask=mask --ld=stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=${flameoMethod[$k]} &> pooling.txt
			else
				$FSLDIR/bin/flameo --cope=cope --vc=varcope --mask=mask --ld=stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=${flameoMethod[$k]} &> pooling.txt
			fi
		# Go back to temp folder
		cd "${temp}"
		done


	# Clean up temp folder
	cd "${temp}"
	rm *
	# Re-initialize variables
	unset -v 'ConFile'
	unset -v 'VarConFile'
	unset -v 'subjects'
	# Go back to studies directory
	cd "${Studies}"

done

#S4POOL
fi





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step five: group analysis.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S5GROUPPOOL is TRUE then execute (jumps to fi)
if [ ${5} = TRUE ] ; then


echo STARTING GROUP ANALYSIS

cd "${WD}"
# Subglobal variables
subjects=0
ConFile=
VarConFile=


# Go to "Group Analysis" folder
cd "${GroupAnalysis}"


# Copy the file groupSubjects.txt from the folder Studies to here
scp -r "${Studies}"/groupSubjects.txt .


#-------------------
#---------# Step one
# Start copy the correct subjects from the data directory and pool them using FLAMES 1+2
echo .... Pooling all subjects "for" group analysis"!!!" ....
	# Read in subjects
	IFS=$'\n' read -d '' -r -a subjects < groupSubjects.txt

	# Go to the temp foler located at the WD folder. Copy the data from the subjects to it.
	cd "${temp}"
	NumSub=${#subjects[@]}
	NumSub=$(($NumSub - 1))

	for j in $(eval echo "{0..$NumSub}"); do
		SubID=${subjects[$j]}
		scp -r $data/$SubID/SessionB/EPI_global/swea/con_0023_$SubID.nii .
		scp -r $data/$SubID/SessionB/EPI_global/swea/varcon_0023_$SubID.nii .

		# Copy the universal mask to the temp folder, rename it with a subject index
		scp -r $data/AllSubMask.nii .
		mv AllSubMask.nii AllSubMask_$SubID.nii

		ConFile=("${ConFile[@]}" con_0023_$SubID.nii)
		VarConFile=("${VarConFile[@]}" varcon_0023_$SubID.nii)
	done

	# Use fslmerge with output again to group analysis folder and then extract the files.
	$FSLDIR/bin/fslmerge -t "$GroupAnalysis/cope" ${ConFile[@]}
	$FSLDIR/bin/fslmerge -t "$GroupAnalysis/varcope" ${VarConFile[@]}
	cd "${GroupAnalysis}"
		gunzip cope.nii.gz
		gunzip varcope.nii.gz

	# Move mask to here
	scp -r "${temp}"/AllSubMask_* .

	# Now create the needed files for flameo (design.mat, design.con and design.grp and a mask) through R file designCrossVal.R
		# Also create the same amount of masks as there are subjects (these have to be pooled using fslmerge before flameo)
		Rscript "${SCRPT}"/designCrossVal.R "$GroupAnalysis/" "$(($NumSub + 1))" &> ROutput_designCrossVal.txt

	# Merge the masks: 4D mask
	$FSLDIR/bin/fslmerge -t "${GroupAnalysis}"/mask.nii AllSubMask_*.nii

	# Remove individual masks
	rm -r AllSubMask_*.nii

	# Now perform flameo to pool subjects: cope is a 4D data structure (4th dimension are the subjects)
	$FSLDIR/bin/flameo --cope=cope --vc=varcope --mask=mask --ld=stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame12&> pooling.txt

	# Unzip cope

# Clean up temp folder
cd "${temp}"
rm *
# Re-initialize variables
unset -v 'ConFile'
unset -v 'VarConFile'
unset -v 'subjects'


# S5GROUPPOOL
fi


#-------------------
#---------# Step two: inference


# If S5GROUPPOOL is TRUE then execute (jumps to fi)
if [ ${6} = TRUE ] ; then


# Cluster based analysis (uncorrected)
if [ "$threshGroup" = cluster ] ; then
	echo .... Cluster wise inference "for" group analysis, uncorrected"!!!" ....
	cd "${GroupAnalysis}"

		echo $($FSLDIR/bin/fslnvols cope) - 1 | bc -l  > stats/dof
		/bin/rm -f stats/zem* stats/zols* stats/mask* # zem and zols are unknown actually... [taken from Nichols blog, so I leave it here]
		$FSLDIR/bin/smoothest -d $(cat stats/dof) -m mask -r stats/res4d > stats/smoothness
		#rm -f stats/res4d*								# Why remove these? Only remove them if you are sure program works and pooling subjects is done (otherwise we need to rerun that and it takes long to do it)!!
		awk '/VOLUME/ {print $2}' stats/smoothness > thresh_zstat1.vol
		awk '/DLH/ {print $2}' stats/smoothness > thresh_zstat1.dlh
		$FSLDIR/bin/fslmaths stats/zstat1 -mas mask thresh_zstat1
		$FSLDIR/bin/cluster -i thresh_zstat1 -c stats/cope1 -t 2.5 -p 0.001 -d $(cat thresh_zstat1.dlh) --volume=$(cat thresh_zstat1.vol) --othresh=thresh_zstat1 -o cluster_mask_zstat1 --connectivity=26 --minclustersize=10 --mm --olmax=lmax_zstat1_STD.txt > cluster_zstat1_std.txt


	# unzip thresh_zstat1 and remove the .nii.gz
	gunzip thresh_zstat1.nii.gz
	# Remove unzipped copes, if they are present (flameo doesn't want two sets of copes in the same folder)
	if [ -f "thresh_zstat1.nii.gz" ]
		then
		rm -r thresh_zstat1.nii.gz
	fi
# Cluster fi
fi



# Voxelwise analysis: FDR corrected.
if [ "$threshGroup" = FDR ] ; then
	echo .... Voxel wise inference "for" group analysis, corrected "for" multiple testing"!!!" ....
	cd "${GroupAnalysis}"

	# Create p-value from the Zstat image.
	$FSLDIR/bin/fslmaths stats/zstat1.nii.gz -ztop stats/pval

	# Now use this pval image to correct for MT
	Rscript "${SCRPT}"/voxelInference.R "$GroupAnalysis" "$GroupAnalysis/stats" "$signLevel" &> ROutput_voxelInference.txt

	# Clean up main and stats folder
	rm cope.nii
	rm varcope.nii
	rm stats/weights1.nii.gz
	rm stats/cope1.nii.gz
	rm stats/mean_random_effects_var1.nii.gz
	rm stats/pe1.nii.gz
	rm stats/p1.nii.gz
	rm stats/res4d.nii.gz
	rm stats/tdof_t1.nii.gz
	rm stats/varcope1.nii.gz
	rm stats/zflame1lowertstat1.nii.gz
	rm stats/zflame1uppertstat1.nii.gz
	rm stats/mask.nii.gz

# FDR fi
fi




# Voxelwise FSL FDR corrections: SAME RESULT AS IN R!
if [ "$threshGroup" = FSLFDR ] ; then
	echo .... Voxel wise inference "for" group analysis, corrected "for" multiple testing"!!!" ....
	cd "${GroupAnalysis}"

	cd stats

	# Create p-value from the Zstat image.
	$FSLDIR/bin/fslmaths zstat1.nii -ztop p1

	# Computes the FDR threshold
	$FSLDIR/bin/fdr -i p1 -m ../mask -q "$signLevel" > FDRPtresh.txt
	Rscript "${SCRPT}"/ptoZ.R "$GroupAnalysis/stats"

	cd ..

	# Get degrees of freedom and smoothness
	echo $($FSLDIR/bin/fslnvols cope) - 1 | bc -l  > stats/dof
	$FSLDIR/bin/smoothest -d $(cat stats/dof) -m mask -r stats/res4d > stats/smoothness

	# Read Zthresh.txt file in which the Z threshold is given
	FDRZthresh=$(echo $cat | awk 'NR==2' stats/FDRZtresh.txt)

	# Make clusters for local maxima without correcting for testing (since we use FDR)
	$FSLDIR/bin/cluster --zstat=stats/zstat1.nii --zthresh="$FDRZthresh" --othresh=thresh_zstat1 > cluster_zstat.txt
	gunzip thresh_zstat1.nii.gz

# FSLFDR fi
fi


#-------------------
#---------# Step three: headers

# Copy header information from subject 33 (random subject): same header, assume safe to copy.
scp -r $data/33/SessionB/EPI_global/swea/con_0023_33.nii .
$FSLDIR/bin/fslcpgeom con_0023_33.nii thresh_zstat1.nii.gz
# Also for tstat1.nii, because we use it later on in transformations
		cd stats
		scp -r $data/33/SessionB/EPI_global/swea/con_0023_33.nii .
		$FSLDIR/bin/fslcpgeom con_0023_33.nii tstat1.nii.gz
		rm -r con_0023_33.nii
		cd ..
rm -r con_0023_33.nii



#S5GROUPINF
fi






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step six: inference in individual studies
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# If S6INFSTUD is TRUE then execute (jumps to fi)
if [ ${7} = TRUE ] ; then

# Go back to studies directory
cd "${Studies}"

# For loop over all studies and then each method of pooling the subjects
for i in $(eval echo "{1..$NS}"); do
	echo Thresholding study:....... $i
	cd Study_$i/
	WRITE=$Studies/Study_$i/

	# For loop over all methods
	for k in $(eval echo "{0..$NumPoolMet}"); do
		MET=${PoolMeth[$k]}
		echo On to method $MET
		cd "${WRITE}${MET}"

		# Cluster based analysis (uncorrected)
		if [ "$threshStudy" = cluster ] ; then

			gunzip mask.nii.gz
			# Remove unzipped masks, if they are present
			if [ -f "mask.nii.gz" ]
				then
        		rm -r mask.nii.gz
			fi

			# Remove unzipped zstats, if they are present
			if [ -f "thresh_zstat1.nii.gz" ]
				then
        		rm -r thresh_zstat1.nii.gz
			fi

			echo $($FSLDIR/bin/fslnvols cope) - 1 | bc -l  > stats/dof
			/bin/rm -f stats/zem* stats/zols* stats/mask* # zem and zols are unknown actually... [taken from Nichols blog, so I leave it here]
			$FSLDIR/bin/smoothest -d $(cat stats/dof) -m mask -r stats/res4d > stats/smoothness
			#rm -f stats/res4d*								# Why remove these? Only remove them if you are sure program works and pooling subjects is done!!
			awk '/VOLUME/ {print $2}' stats/smoothness > thresh_zstat1.vol
			awk '/DLH/ {print $2}' stats/smoothness > thresh_zstat1.dlh
			$FSLDIR/bin/fslmaths stats/zstat1 -mas mask thresh_zstat1
			$FSLDIR/bin/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -p $signLevelStud -d $(cat thresh_zstat1.dlh) --volume=$(cat thresh_zstat1.vol) --othresh=thresh_zstat1 -o cluster_mask_zstat1 --connectivity=26 --minclustersize=$minclustersizeStud --mm --olmax=lmax_zstat1_STD.txt > cluster_zstat1_std.txt


			gunzip thresh_zstat1.nii.gz

		# Cluster fi
		fi


		# Peak level analysis (FDR corrected). But report local maxima within cluster.
		if [ "$threshStudy" = peak ] ; then

			echo $($FSLDIR/bin/fslnvols cope) - 1 | bc -l  > stats/dof
			/bin/rm -f stats/zem* stats/zols* stats/mask* # zem and zols are unknown actually... [taken from Nichols blog, so I leave it here]
			$FSLDIR/bin/smoothest -d $(cat stats/dof) -m mask -r stats/res4d > stats/smoothness
			awk '/VOLUME/ {print $2}' stats/smoothness > thresh_zstat1.vol
			awk '/DLH/ {print $2}' stats/smoothness > thresh_zstat1.dlh
			cd stats
			gunzip zstat1.nii.gz
			cd ..
			gunzip mask.nii.gz


			# R script for finding adjusted Z-threshold
			Rscript "${SCRPT}"/peakInference.R "$WRITE/$MET" "$WRITE/$MET/stats" "$excursionSet" "$signLevelStud" &> ROutput_peakInference.txt
			# Read Zthresh.txt file in which the Z threshold is given
			read -d '' -r -a Zthreshold < Zthresh.txt
			$FSLDIR/bin/cluster --zstat=stats/zstat1.nii --zthresh="$Zthreshold" --othresh=thresh_zstat1 > cluster_zstat.txt
			# Unzip thresholded map, if not unzipped
			if [ -f "thresh_zstat1.nii.gz" ]
				then
        		gunzip thresh_zstat1.nii.gz
			fi

			# Make a text file: lmax_zstat1_STD.txt which contains columns: Cluster Index, Value, x, y, z from clusters with minimum cluster size = minclustersizeStud
					# So we do not look at clusters with cluster size < minclustersizeStud
			Rscript "${SCRPT}"/lmax_zstat.R "$WRITE/$MET" "$minclustersizeStud"

		# Peak fi
		fi


		# FDR correction: Nichols. Report local maxima of SPM.
		if [ "$threshStudy" = FDR ] ; then

			cd stats
			# Create p-value from the Zstat image.
			$FSLDIR/bin/fslmaths zstat1.nii -ztop p1

			# Computes the FDR threshold
			$FSLDIR/bin/fdr -i p1 -m ../mask -q 0.05 > FDRPtresh.txt
			Rscript "${SCRPT}"/ptoZ.R "$WRITE/$MET/stats" &> ROutput_tmp.txt
			rm ROutput_tmp.txt

			cd ..

			# Read Zthresh.txt file in which the Z threshold is given
			FDRZthresh=$(echo $cat | awk 'NR==2' stats/FDRZtresh.txt)

			# Make clusters for local maxima without correcting for testing (since we use FDR)
			$FSLDIR/bin/cluster --zstat=stats/zstat1.nii --zthresh="$FDRZthresh" --othresh=thresh_zstat1 > cluster_zstat.txt

			# Unzip thresholded map, if not unzipped
			if [ -f "thresh_zstat1.nii.gz" ]
				then
        		gunzip thresh_zstat1.nii.gz
			fi

			# Make a text file: lmax_zstat1_STD.txt which contains columns: Cluster Index, Value, x, y, z from clusters with minimum cluster size = minclustersizeStud
					# So we do not look at clusters with cluster size < minclustersizeStud
			Rscript "${SCRPT}"/lmax_zstat.R "$WRITE/$MET" "$minclustersizeStud"

			# Clean up main and stats folder
			rm cope.nii
			rm varcope.nii
			rm stats/weights1.nii.gz
			#rm stats/cope1.nii.gz
			rm stats/mean_random_effects_var1.nii.gz
			rm stats/pe1.nii.gz
			rm stats/p1.nii.gz
			rm stats/res4d.nii.gz
			rm stats/tdof_t1.nii.gz
			#rm stats/varcope1.nii.gz
			rm stats/zflame1lowertstat1.nii.gz
			rm stats/zflame1uppertstat1.nii.gz
			rm stats/mask.nii.gz

		# Peak fi
		fi


	# Copy header information from subject 33 (random subject): same header, assume safe to copy.
	scp -r $data/33/SessionB/EPI_global/swea/con_0023_33.nii .
	$FSLDIR/bin/fslcpgeom con_0023_33.nii thresh_zstat1.nii.gz

	rm -r con_0023_33.nii


	done
	cd "${Studies}"

done
cd "${WD}"


# S6INFSTUD
fi




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step seven: convert all (group if cluster analysis and individual studies) foci to MNI
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# If S7FOCICONV is TRUE then execute (jumps to fi)
if [ ${8} = TRUE ] ; then


# If cluster analysis is used for group study: convert foci locations
if [ "$threshGroup" = cluster ] ; then
	echo 'Need to code this part -- SKIPPING convertions of foci'
fi

if [ "$threshGroup" = FDR ] ; then
	echo 'No conversion needed for group analysis -- going to individual studies'
fi

# Go to directory
cd "${Studies}"

# For loop over all studies and then each method of pooling the subjects
for i in $(eval echo "{1..$NS}"); do
	echo Converting study:....... $i
	cd Study_$i/
	WRITE=$Studies/Study_$i/

	# For loop over all methods
	for k in $(eval echo "{0..$NumPoolMet}"); do
		MET=${PoolMeth[$k]}
		cd "${WRITE}${MET}"

		# Starting R script convertFoci.R
		Rscript "${SCRPT}"/convertFoci.R "$WRITE/$MET" "lmax_zstat1_STD.txt" &> ROutput_convertFoci.txt

	done
	cd "${Studies}"

done
# Back to WD
cd "${WD}"



# S7FOCICONV
fi




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step eight: perform meta-analyses
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S8METAN is TRUE then execute (jumps to fi)
if [ ${9} = TRUE ] ; then


cd "${MetaAnalyses}"

# SEED has been reset in obtaining the results while coding parts of the program, hence for reproducibility, we will leave it here.
SEED=$11121990

# Subglobal variables
	# Array of used methods:
	MetaMeth=(	FixedEffUn
			 	RanEffUn)
	MetaMethod=(Fixed
				Random)


	# Number of meta-analyses without the pooling
	NumMeta=${#MetaMeth[@]}
	NumMeta=$(($NumMeta - 1))

# Create the folders
for i in $(eval echo "{0..$NumMeta}"); do
	mkdir "${MetaMeth[$i]}"
	cd "${MetaMeth[$i]}"
	for k in $(eval echo "{0..$NumPoolMet}"); do
		mkdir "${PoolMeth[$k]}"
	done
	cd ..
done
cd "${MetaAnalyses}"


# For loop over all methods of the meta-analyses,
for i in $(eval echo "{0..$NumMeta}"); do
	cd "${MetaMeth[$i]}"
	# and then the pooling methods
	for k in $(eval echo "{0..$NumPoolMet}"); do
		cd "${PoolMeth[$k]}"
		CurrentMethod=("$MetaAnalyses"/"${MetaMeth[$i]}"/"${PoolMeth[$k]}"/)

		# Take the peak (local maxima) locations and put them in a file with the right format (number of studies for loop).
		echo .... Writing "local" maxima 'for' all studies .....
		for j in $(eval echo "{1..$NS}"); do
			LocStud=("${Studies}"/Study_$j/"${PoolMeth[$k]}"/)
			NameStud=(MNI_lmax_zstat1_STD.txt)

				# CAREFUL NOT TO MESS WITH THE ORDER OR THE ARGUMENTS (own argument is to denote not to write to ALE format but own format)
			Rscript "${SCRPT}"/preparePeaks.R "$WD" "$j" "${arrSub[$(($j - 1))]}" "$LocStud" "$NameStud" "$CurrentMethod" "peakFileStudy_" "Own" "FALSE" "6"	&> ROutput_preparePeaks.txt
		done

		# Perform the meta-analysis and threshold based on ThrSetMeth variable
			echo .... Performing Meta-Analysis .... Calculate P-Values ....
		Rscript "${SCRPT}"/metaAnalysis.R "$WD" "$SEED" "${MetaMethod[$i]}" "$locMask" "$NS" "$CurrentMethod" "peakFileStudy_" "$CurrentMethod" "$COMPUTER" &> ROutput_metaAnalysis.txt


		cd ..
	done
		cd "${MetaAnalyses}"
done



# S8METAN
fi




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step nine: ALE meta-analysis.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S9ALE is TRUE then execute (jumps to fi)
if [ ${10} = TRUE ] ; then

echo "Starting preparation of ALE meta-analyses"

# Subglobal variables
	# Array of used methods:
	ALEMetaMeth=(ALEUn
				ALECluster
				ALEpN
			 	ALEpID)


	# Number of ALE meta-analyses (not considering the pooling methods)
	NumALEMeta=${#ALEMetaMeth[@]}
	NumALEMeta=$(($NumALEMeta - 1))


# Go to meta analysis folder
cd "${MetaAnalyses}"

# Make ALE folder
mkdir ALE
cd ALE

# For loop  over pooling methods
for k in $(eval echo "{0..$NumPoolMet}"); do
	# folder according to pooling method
	mkdir "${PoolMeth[$k]}"
	cd "${PoolMeth[$k]}"
	# DataRaw
	mkdir DataRaw

	# %%%% Now run the preparePeaks function %%%%
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	echo %%%%%%%%% Writing Local Maxima %%%%%%%%%
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	# Where to write to:
	OutputFociALE=("$MetaAnalyses"/ALE/"${PoolMeth[$k]}"/DataRaw/)

	# Take the peak (local maxima) locations and put them in a file with the right format (number of studies for loop).
	for j in $(eval echo "{1..$NS}"); do
		LocStud=("${Studies}"/Study_$j/"${PoolMeth[$k]}"/)
		NameStud=(MNI_lmax_zstat1_STD.txt)

			# CAREFUL NOT TO MESS WITH THE ORDER OF THE ARGUMENTS
		Rscript "${SCRPT}"/preparePeaks.R "$WD" "$j" "${arrSub[$(($j - 1))]}" "$LocStud" "$NameStud" "$OutputFociALE" "peakFileStudy_" "ALE" "FALSE" "6" &> ROutput_preparePeaks.txt
	done

	# Go back to ALE folder
	cd "${MetaAnalyses}"/ALE

done

# S9ALE
fi





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step ten: transformations of ALE data
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# If S10TRANALE is TRUE then execute (jumps to fi)
if [ ${11} = TRUE ] ; then

if [ "$MATLABVERSION" = NO ] ; then
	echo "Unfortunately, the MATLAB version of ALE does not work on the HPC infrastructure.
				Hence you need to switch to local PC/MAC/server and run it there with a
				separate script."

# fi for MATLAB
fi

# S10TRANALE
fi




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step eleven: group study of subjects in the MA
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# If S11GROUPMA is TRUE then execute (jumps to fi)
if [ ${12} = TRUE ] ; then


cd "${WD}"
echo STARTING A GROUP ANALYSIS OF ALL SUBJECTS USED IN THE META-ANALYSIS: ARROW A

# Subglobal variables
subjects=0
ConFile=
VarConFile=

#-------------------
#---------# Step one: get all data in temp folder, have an array (ConFile and VarConFile with all subjects) and pool the subjects

# For loop over all studies
cd "${Studies}"
for i in $(eval echo "{1..$NS}"); do
	cd Study_$i/
	WRITE=$SubMetGroupAn
	ConFile=
	VarConFile=

	# Get subject ID's
	IFS=$'\n' read -d '' -r -a subjects < study_$i.txt
	# ID variable ONLY for looping purpose (if number of subjects = 12: then BASH needs to loop from 0-11 !!!)
	NumSub=${#subjects[@]}
	NumSub=$(($NumSub - 1))

	# Go to the temp foler located at the data folder. Copy the data from the subjects to it.
	cd "${temp}"
	for j in $(eval echo "{0..$NumSub}"); do
	SubID=${subjects[$j]}
	scp -r $data/$SubID/SessionB/EPI_global/swea/con_0023_$SubID.nii .
	scp -r $data/$SubID/SessionB/EPI_global/swea/varcon_0023_$SubID.nii .

	# Copy the universal mask to the temp folder, rename it with a subject index
	scp -r $data/AllSubMask.nii .
	mv AllSubMask.nii AllSubMask_$SubID.nii

	ConFile=("${ConFile[@]}" con_0023_$SubID.nii)
	VarConFile=("${VarConFile[@]}" varcon_0023_$SubID.nii)
	done

	# Back to Studies
	cd "${Studies}"

	# Reset subjects (NOT ConFile and VarConFile!)
	unset -v 'subjects'
done


# Then we pool all subjects in the same fashion of the GT
	cd "${temp}"

	# Use fslmerge with output again to group analysis folder and then extract the files.
	$FSLDIR/bin/fslmerge -t "$SubMetGroupAn/cope" ${ConFile[@]}
	$FSLDIR/bin/fslmerge -t "$SubMetGroupAn/varcope" ${VarConFile[@]}
	cd "${SubMetGroupAn}"
		gunzip cope.nii.gz
		gunzip varcope.nii.gz

	# Move mask to here
	scp -r "${temp}"/AllSubMask_* .

	# Now create the needed files for flameo (design.mat, design.con and design.grp and a mask) through R file designCrossVal.R
		# Also create the same amount of masks as there are subjects (these have to be pooled using fslmerge before flameo)
		Rscript "${SCRPT}"/designCrossVal.R "$SubMetGroupAn/" "$groupNS" &> ROutput_designCrossValSubMetGroupAn.txt

	# Merge the masks: 4D mask
	$FSLDIR/bin/fslmerge -t "${SubMetGroupAn}"/mask.nii AllSubMask_*.nii

	# Remove individual masks
	rm -r AllSubMask_*.nii

	echo 'Pooling the subjects'
	# Now perform flameo to pool subjects: cope is a 4D data structure (4th dimension are the subjects)
	$FSLDIR/bin/flameo --cope=cope --vc=varcope --mask=mask --ld=stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame12&> poolingSubMetGroupAn.txt

	# Unzip cope

# Clean up temp folder
cd "${temp}"
rm *
# Re-initialize variables
unset -v 'ConFile'
unset -v 'VarConFile'
unset -v 'subjects'




#-------------------
#---------# Step two: inference
echo 'Inference'
# Voxelwise analysis: FDR corrected.
if [ "$threshGroup" = FDR ] ; then
	echo .... Voxel wise inference "for" group analysis, corrected "for" multiple testing"!!!" ....
	cd "${SubMetGroupAn}"

	# Create p-value from the Zstat image.
	$FSLDIR/bin/fslmaths stats/zstat1.nii.gz -ztop stats/pval

	# Now use this pval image to correct for MT
	Rscript "${SCRPT}"/voxelInference.R "$SubMetGroupAn" "$SubMetGroupAn/stats" "$signLevel" &> ROutput_voxelInferenceSubMetGroupAn.txt

	# Clean up main and stats folder
	rm cope.nii
	rm varcope.nii
	rm stats/weights1.nii.gz
	rm stats/cope1.nii.gz
	rm stats/mean_random_effects_var1.nii.gz
	rm stats/pe1.nii.gz
	rm stats/p1.nii.gz
	rm stats/res4d.nii.gz
	rm stats/tdof_t1.nii.gz
	rm stats/varcope1.nii.gz
	rm stats/zflame1lowertstat1.nii.gz
	rm stats/zflame1uppertstat1.nii.gz
	rm stats/mask.nii.gz

# FDR fi
fi


#-------------------
#---------# Step three: headers

# Copy header information from subject 33 (random subject): same header, assume safe to copy.
scp -r $data/33/SessionB/EPI_global/swea/con_0023_33.nii .
$FSLDIR/bin/fslcpgeom con_0023_33.nii thresh_zstat1.nii.gz
# Also for tstat1.nii
		cd stats
		scp -r $data/33/SessionB/EPI_global/swea/con_0023_33.nii .
		$FSLDIR/bin/fslcpgeom con_0023_33.nii tstat1.nii.gz
		rm -r con_0023_33.nii
		cd ..
rm -r con_0023_33.nii


# S11GROUPMA
fi



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Appendix: remove the temp folder
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Go to RUNWE and remove the folder
cd "${RUNWD}"
rm -r temp
