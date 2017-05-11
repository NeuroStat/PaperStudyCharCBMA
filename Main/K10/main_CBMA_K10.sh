#!/bin/sh

####################
#### TITLE:     -------------K = 10 VERSION------------
####			Code to sample 2 x 400 subjects over all studies from the Imagen data pool
####			TEST condition: subsample batch of 400 into studies and perform fixed effects, ols or mixed effects group analyses.
####			And then perform the ALE, fixed and random effects meta-analysis.
####			EVALUATION condition: batch of 400 into one mixed effects group analysis.
####			This is repeated for several unique iterations.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_CBMA/PaperStudyCharCBMA.git/Main/K10
#### First Modified: 13/11/2014
#################


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# INFORMATION
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Structure of folders after running this:
# Scripts (location of the main program and all scripts it uses)
#--# Working Directory (called Results were studies (after pooling), meta analyses and group analysis are stored in each time a folder for a run)
#----# The number of RUN we are in (e.g. Run_1)
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

# Actual IMAGEN (or other input) data can be in any location!


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DIRECTORIES
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Location of the data
data="<<LocationOfData>>"

# Location of the ANONYMIZED identifiers for subject and scanner location
AllIDs=$data/OurIDs

# Location of scripts folder
SCRPT=PaperStudyCharCBMA/Main/K10


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SELECT WHICH PART OF SCRIPT TO EXECUTE
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S2PREPARATION=TRUE
S3SAMPLING=TRUE
S4POOL=TRUE
S5GROUPPOOL=TRUE
S6GROUPINF=TRUE
S7INFSTUD=TRUE
S8FOCICONV=TRUE
S9METAN=TRUE





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step one: Global variables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Empty array of amount of subjects in each study
arrSub=()
# Location of mask for the meta-analysis
locMask=($data/AllSubMask.nii)
# Number of studies in the meta-analysis
NS=10
# Number of subjects in the group analysis
groupNS=200
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
# In which iteration/run are we?
# We have run this on high performance cluster which takes 1-7 as argument.
# Hence this is run with value 1, 2, 3, 4, 5, 6, 7
RUN=1



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step two: preparation
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S2PREPARATION is TRUE then execute
if [ $S2PREPARATION = TRUE ] ; then

# Create Results folder
cd "${SCRPT}"
# If we are in the first run, then create folder results
if [ "$RUN" = 1 ] ; then
	mkdir Results
fi
RUNWD=$SCRPT/Results

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
# Make the directories corresponding to each study, then go into the folder and make
# 			fixed effects, OLS, mixed effects pooling

cd "${Studies}"
for i in $(eval echo "{1..$NS}"); do
 mkdir Study_$i
 cd Study_$i
 mkdir FixedEffect
 mkdir OLS
 mkdir MixedEffect
 cd ..
done

# end of S2PREPARATION
fi



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step three: Sampling
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S3SAMPLING is TRUE then execute (jumps to fi)
if [ $S3SAMPLING = TRUE ] ; then

echo Sampling subjects
cd "${SCRPT}"
Rscript "${SCRPT}"/LoadSubjK10.R "$Studies" "$RUN" "$SCRPT" &> ROutput_SamplingCrossVal.txt


# S3SAMPLING
fi





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step four: pooling subjects in studies
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Always execute this part of the code!

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# If S4POOL is TRUE then execute
if [ $S4POOL = TRUE ] ; then

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
		scp -r $data/$SubID/con_0023_$SubID.nii .
		scp -r $data/$SubID/varcon_0023_$SubID.nii .
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
# Step five: pooling subjects in evaluation condition (reference)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S5GROUPPOOL is TRUE then execute (jumps to fi)
if [ $S5GROUPPOOL = TRUE ] ; then


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
		scp -r $data/$SubID/con_0023_$SubID.nii .
		scp -r $data/$SubID/varcon_0023_$SubID.nii .

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




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step six: inference in evaluation condition (reference)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S6GROUPINF is TRUE then execute (jumps to fi)
if [ $S6GROUPINF = TRUE ] ; then

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

#-------------------
#---------# Step 2: need to copy headers

# Copy header information from subject 33 (random subject): same header, assume safe to copy.
scp -r $data/33/con_0023_33.nii .
$FSLDIR/bin/fslcpgeom con_0023_33.nii thresh_zstat1.nii.gz
# Also for tstat1.nii, because we use it later on in transformations
		cd stats
		scp -r $data/33/con_0023_33.nii .
		$FSLDIR/bin/fslcpgeom con_0023_33.nii tstat1.nii.gz
		rm -r con_0023_33.nii
		cd ..
rm -r con_0023_33.nii



#S6GROUPINF
fi





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step seven: inference in individual studies
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S7INFSTUD is TRUE then execute (jumps to fi)
if [ $S7INFSTUD = TRUE ] ; then

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

		# FDR correction. Report local maxima of SPM.
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
			rm stats/mean_random_effects_var1.nii.gz
			rm stats/pe1.nii.gz
			rm stats/p1.nii.gz
			rm stats/res4d.nii.gz
			rm stats/tdof_t1.nii.gz
			rm stats/zflame1lowertstat1.nii.gz
			rm stats/zflame1uppertstat1.nii.gz
			rm stats/mask.nii.gz

		# FDR fi
		fi

	# Copy header information from random subject: same header, assume safe to copy.
	scp -r $data/33/con_0023_33.nii .
	$FSLDIR/bin/fslcpgeom con_0023_33.nii thresh_zstat1.nii.gz

	rm -r con_0023_33.nii

	done
	cd "${Studies}"

done
cd "${WD}"


# S7INFSTUD
fi





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step eight: convert foci from individual studies to MNI
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# We use a transformation from the processed data to MNI space.
# See convertFoci.R

# If S8FOCICONV is TRUE then execute (jumps to fi)
if [ $S8FOCICONV = TRUE ] ; then

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

# S8FOCICONV
fi



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step nine: perform effect size based meta-analyses
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# If S9METAN is TRUE then execute (jumps to fi)
if [ $S9METAN = TRUE ] ; then

cd "${MetaAnalyses}"

# SEED has been reset in obtaining the results while coding parts of the program, hence for reproducibility, we will leave it here.
# New data with same seed, no problem.
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
		Rscript "${SCRPT}"/metaAnalysis.R "$WD" "$SEED" "${MetaMethod[$i]}" "$locMask" "$NS" "$CurrentMethod" "peakFileStudy_" "$CurrentMethod" &> ROutput_metaAnalysis.txt

		cd ..
	done
		cd "${MetaAnalyses}"
done

# S9METAN
fi


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step eleven: ALE meta-analysis.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Seperate scripts have been run to execute the ALE meta-analyses.
# These scripts were provided by prof. dr. Simon Eickhoff (personal communication).
# We have not explicitly asked for permission to include these in a public repository.
# Hence if interest, please contact Simon.



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Appendix: remove the temp folder
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Go to RUNWE and remove the folder
cd "${RUNWD}"
rm -r temp













