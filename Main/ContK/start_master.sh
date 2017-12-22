#!/bin/sh
#
#
#PBS -N FixRanContKFast
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=11:00:00
#PBS -l vmem=30GB
#


#----------------------------------------------------#
# MODULES TO LOAD IN
module load R/3.2.3-intel-2016a
module load FSL/5.0.9-intel-2016a
module swap Java Java/1.8.0_92
. $FSLDIR/etc/fslconf/fsl.sh
#----------------------------------------------------#


#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#

#----------------------------------------------------#
# PARAMETERS FOR IN MASTER FILE: WHICH STEPS TO RUN?
S2FOLDERS=FALSE
S3SAMPLING=TRUE
S4POOL=FALSE
S5GROUPPOOL=FALSE
S5GROUPINF=FALSE
S6INFSTUD=FALSE
S7FOCICONV=FALSE
S8METAN=TRUE
S9ALE=FALSE
S10TRANALE=FALSE		# DEFAULT FALSE: is obsolete now
S11GROUPMA=FALSE		# DEFAULT FALSE: is obsolete now
#----------------------------------------------------#

#----------------------------------------------------#
# LOCATION OF SCRIPT TO RUN
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/ContKStudiesFast
cd $srcdir
#----------------------------------------------------#

#----------------------------------------------------#
# SELECT NUMBER OF K
KSEL=(12 12 12 12 12 14 14 14 14 14 16 16 16 16 18 18 18 30 30)
K=${KSEL[${PBS_ARRAYID}]}
#----------------------------------------------------#

#----------------------------------------------------#
# SELECT NUMBER OF ASSOCIATED RUNS (FOLDS). See file SamplingCVContK.R for numbers.
RUNSEL=(1 2 3 4 5 1 2 3 4 5 1 2 3 4 1 2 3 1 2)
RUN=${RUNSEL[${PBS_ARRAYID}]}
#----------------------------------------------------#

#----------------------------------------------------#
# SELECT NUMBER OF SUBJECTS IN GT. See file SamplingCVContK.R for numbers.
NTEST=(240 240 240 240 240 280 280 280 280 280 320 320 320 320 360 360 360 600 600)
N=${NTEST[${PBS_ARRAYID}]}
#----------------------------------------------------#

#----------------------------------------------------#
# GO TIME: CAREFUL WITH ORDER OF ARGUMENTS!
	./master_fixran_ContK.sh "$RUN" "$S2FOLDERS" "$S3SAMPLING" "$S4POOL" "$S5GROUPPOOL" "$S5GROUPINF" "$S6INFSTUD" "$S7FOCICONV" "$S8METAN" "$S9ALE" "$S10TRANALE" "$S11GROUPMA" "$vsc" "$K" "$N"
#----------------------------------------------------#




echo "job finished"
