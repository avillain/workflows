#!/bin/bash

#$ -S /bin/bash

#$ -cwd

#$ -q pf4 

#$ -N freec_test

#$ -l mem_total=25G

module load samtools
module load freec
module load R

####################
# PRECONFIGURATION #
####################

#### $1 : samfile.sam $2 : reference.bam $3 : reference.len

#directories
TMP="./tmp"
RESULTS="./results"

if [[ ! -s $RESULTS ]]
     then
          mkdir $RESULTS
fi

if [[ ! -s $TMP ]]
     then
          mkdir $TMP
fi

#filenames
prefix=`basename ${1%.*}`
lenfile=$3
bamfile=$TMP"/"$prefix".bam"
bamctrl=$2
config=$TMP"/"$prefix".config"

echo "[general]
chrLenFile = $lenfile
minExpectedGC=0.25
minExpectedGC=0.45
window = 3000

step = 1000
ploidy = 1
#breakPointThreshold = -.001

telocentromeric = 10000

GCcontentProfile = GC_profile.cnp

intercept=1
minMappabilityPerWindow = 0.7

outputDir = $RESULTS 

breakPointType=4

[sample]

mateFile = $bamfile

inputFormat = BAM
mateOrientation = FR

[control]

mateFile = $bamctrl

inputFormat = BAM
mateOrientation = FR" > $config


#samtobam
if [ ! -f $bamfile ]
then
    samtools view -bS $1 > $bamfile
fi


#freec
freec -conf $config

cat /pasteur/homes/avillain/scripts/processing/makeGraph_spombe1.sh | R --slave --args 1 $RESULTS"/"`basename $bamfile`"_I_II_III_ratio.txt"
cat /pasteur/homes/avillain/scripts/processing/makeGraph_spombe2.sh | R --slave --args 1 $RESULTS"/"`basename $bamfile`"_MT_MTR_AB_ratio.txt"

