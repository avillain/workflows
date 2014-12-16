#!/bin/bash

#$ -S /bin/bash

#$ -cwd

#$ -q pf4 

#$ -N pindel_test

#$ -l mem_total=25G

module load samtools
module load pindel
module load SHRiMP bwa
module load PRISM
module load SOAPindel
module load tabix
module load snpEff/3.5

####################
# PRECONFIGURATION #
####################


#### $1 : samfile.sam $2 : reference.fa $3 : compare.vcf.gz $4 : snpeff.config $5 : comparepindel.vcf.gz

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
bamfile=$TMP"/"$prefix".bam"
sortedbam=$TMP"/"$prefix"_sorted"
sortedbamfile=$sortedbam.bam
bamindex=$sortedbamfile.bai
reference=$TMP"/"`basename $2`
refindex=$reference.fai
compare=$RESULTS"/"`basename $3`

compindel=$5
cp $2 $reference
cp $3 $compare

#insert size
insertfile=$TMP"/"$prefix"_insertsize.txt"

if [ ! -f $insertfile ]
then
    getinsertsize.py $1 > $insertfile
fi

insertmean=`grep -o -P "span: mean ([0-9]+)" $insertfile | cut -d " " -f 3`
insertsd=`grep -o -P "span: mean [0-9]+\.[0-9]+, STD=([0-9]+)" $insertfile | cut -d "=" -f 2`
readslength=`grep -o -P "length: mean [0-9]+" $insertfile | cut -d " " -f 3`

#samtobam, bam sorting, bam indexing, reference indexing
if [ ! -f $bamfile ]
then
    samtools view -bS $1 > $bamfile
fi

if [ ! -f $sortedbamfile ]
then
    samtools sort $bamfile $sortedbam
fi

if [ ! -f $bamindex ]
then
    samtools index $sortedbamfile $bamindex
fi

if [ ! -f $refindex ]
then
    samtools faidx $reference $refindex
fi


#############
# DETECTION #
#############

### Pindel

pindelconfigfile=$TMP"/""pindel$prefix.config"
pindout=$TMP"/"$prefix"_D"
pindvcf=$TMP"/"$prefix"_pindel.vcf"
pindvcf_filt=$RESULTS"/"$prefix"_pindel_filtered.vcf"
pindelsubtracted=$RESULTS"/"$prefix"_pindel_subtracted.vcf"
pindvcf0=$RESULTS"/"$prefix"_pindel_filt0.vcf"
pindvcf25=$RESULTS"/"$prefix"_pindel_filt25.vcf"
pindvcf50=$RESULTS"/"$prefix"_pindel_filt50.vcf"
pindvcf100=$RESULTS"/"$prefix"_pindel_filt100.vcf"

pindel8cols=$RESULTS"/"$prefix"_pindel_filtered_8cols.vcf"



#configuration
echo -e "$sortedbamfile\t$insertmean\t$prefix" > $pindelconfigfile

#detection
pindel -f $reference -i $pindelconfigfile -c ALL -o $TMP"/"$prefix

#conversion
pindel2vcf -p $pindout -r $reference -R spombe -d 09052011 -v $pindvcf

#filtering
grep "#" $pindvcf > $pindvcf_filt
awk '!/^#.*/ {split($10,tab,","); if(tab[2]>=10) print $0}' $pindvcf >> $pindvcf_filt

grep "^##" $pindvcf_filt > $pindel8cols
awk 'BEGIN {OFS = "\t"} !/^##.*/ { print $1, $2, $3, $4, $5, $6, $7, $8 }' $pindvcf_filt >> $pindel8cols

#subtraction
compindelname=`basename $compindel`
compindelgz=$compindel.gz
pindvcf_filtgz=$pindvcf_filt.gz

if [ ! -f $pindelvcf_filtgz ]
then
    bgzip $pindelvcf_filt
fi

tabix -f -p vcf $pindelvcf_filtgz


if test ${compindel##*.} != 'gz'
then
    bgzip $compindel
    comp2=$compindel.gz
    compindel=$comp2
fi

tabix -f -p vcf $compindel

vcf-isec -f -c $pindelfvcf_filtgz $compindel > $pindelsubtracted

# division in subclasses
grep "#" $pindelsubtracted > $pindvcf0
awk 'BEGIN {OFS = "\t"} !/^#.*/ {split($10,geno,":");split(geno[2],ad,","); if(ad[2]>=10 && ad[2]*100/(ad[1]+ad[2])<=15) print $0}' $pindelsubtracted >> $pindvcf0

grep "#" $pindelsubtracted > $pindvcf25
awk 'BEGIN {OFS = "\t"} !/^#.*/ {split($10,geno,":");split(geno[2],ad,","); if(ad[2]>=10 && ad[2]*100/(ad[1]+ad[2])>15 && ad[2]*100/(ad[1]+ad[2])<=35) print $0}' $pindelsubtracted >> $pindvcf25

grep "#" $pindelsubtracted > $pindvcf50
awk 'BEGIN {OFS = "\t"} !/^#.*/ {split($10,geno,":");split(geno[2],ad,","); if(ad[2]>=10 && ad[2]*100/(ad[1]+ad[2])>35 && ad[2]*100/(ad[1]+ad[2])<=65) print $0}' $pindelsubtracted >> $pindvcf50

grep "#" $pindelsubtracted > $pindvcf100
awk 'BEGIN {OFS = "\t"} !/^#.*/ {split($10,geno,":");split(geno[2],ad,","); if(ad[2]>=10 && ad[2]*100/(ad[1]+ad[2])>65) print $0}' $pindelsubtracted >> $pindvcf100


### SOAPindel

soapvcf=$TMP"/"$prefix"_soap.vcf"
soapvcf_filt=$RESULTS"/"$prefix"_soap_filtered.vcf"
soap8cols=$RESULTS"/"$prefix"_soap_filtered_8cols.vcf"

mappinglist=$TMP"/"$prefix"_mapping.list"
readslist=$TMP"/"$prefix"_reads.list"
searchprefix=$prefix
#searchprefix="1623-1"

#configuration
ls /pasteur/projets/NGS-Dyngen/fastq_files/clean_files/*$searchprefix*.fastq > $readslist
echo -e "$sortedbamfile\t$insertmean\t$insertsd\t$readslength\tPAIR" > $mappinglist

#detection
indel_detection.ibam.pl $mappinglist $reference $readslist -wd $TMP

#merging chromosomes files
grep "#" $TMP/result/I:*/*indel.vcf > $soapvcf
grep -hv "#" $TMP/result/*/*indel.vcf >> $soapvcf

sed -i "s/chrI/I/g" $soapvcf
sed -i "s/chrII/II/g" $soapvcf
sed -i "s/chrMT/MT/g" $soapvcf

#filtering
grep "##" $soapvcf > $soapvcf_filt
grep "#CHR" $soapvcf >> $soapvcf_filt
awk '!/^#.*/ {split($8,tab,";"); split(tab[1],tab2,"="); if(tab2[2]>=5) print $0}' $soapvcf >> $soapvcf_filt

grep "^##" $soapvcf_filt > $soap8cols
awk 'BEGIN {OFS = "\t"} !/^##.*/ { print $1, $2, $3, $4, $5, $6, $7, $8 }' $soapvcf_filt >> $soap8cols

### PRISM (doit être lancé pour chaque chromosome)

prismvcf_filt=$RESULTS"/"$prefix"_prism_filtered.vcf"
refdir=$reference".DB.SPLIT/*"
prism8cols=$RESULTS"/"$prefix"_prism_filtered_8cols.vcf"

#detection
#for i in $refdir; do v=`basename ${i%.*}`; chrname=${v#chr}; j=$TMP"/"$prefix"_"$chrname".sam"; samtools view -h $sortedbamfile $chrname > $j; run_PRISM.sh -m $insertmean -e $insertsd -r $i -i $j -I $TMP -O $TMP"/"$chrname; done
for i in $refdir; do chrname=`basename ${i%.*}`; run_PRISM.sh -m $insertmean -e $insertsd -r $i -i $1 -I $TMP -O $TMP"/"$chrname; done

#filtering
for i in $refdir; do chrname=`basename ${i%.*}`; awk '{ if ($6>=5) print $0 }' $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv" > $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv_supported"; done

#conversion
for i in $refdir; do chrname=`basename ${i%.*}`; out=$TMP"/"$prefix"_"$chrname"_prism.vcf"; prism2vcf.py $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv_supported" $reference $out; done
cat $TMP"/"*prism.vcf > $prismvcf_filt

echo -e "##fileformat=VCFv4.0\n###fileDate=20140724\n###source=prism\n###reference=/pasteur/projets/NGS-Dyngen/large_indels/test/J64-G15/./tmp/pombe_09052011.fasta.DB.SPLIT\n###INFO=<ID=PROGRAM,Number=1,Type=String,Description="Total number of reads in haplotype window">\n###INFO=<ID=SVTYPE,Number=1,Type=String,Description="Total number of reads 2 in haplotype window">\n###INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Left flank length and right flank length">\n###FILTER=<ID=q10,Description="Quality below 10">\n###FILTER=<ID=hp10,Description="Reference homopolymer length was longer than 10">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > $prism8cols
awk 'BEGIN {OFS = "\t"} !/^##.*/ { $2= $2 - 1; print $1, $2, $3, $4, $5, $6, $7, $8 }' $prismvcf_filt >> $prism8cols

###############
# COMPARAISON #
###############

pindgz=$pindel8cols.gz
soapgz=$soap8cols.gz
prismgz=$prism8cols.gz
union=$RESULTS"/"$prefix"_pindel_soap_prism_union.vcf"
compname=`basename $compare .vcf`
minus=$RESULTS"/"$prefix"_pindel_soap_prism_minus_"${compname%.*}".vcf"

bgzip $pindel8cols
bgzip $soap8cols
bgzip $prism8cols
tabix -f -p vcf $pindgz
tabix -f -p vcf $soapgz
tabix -f -p vcf $prismgz

#combinaison des outils
vcf-isec -a -n +1 $pindgz $soapgz $prismgz > $union

if test ${union##*.} != 'gz'
then
    bgzip $union
    un3=$union.gz
    union=$un3
fi

tabix -f -p vcf $union

if test ${compare##*.} != 'gz'
then
    bgzip $compare
    comp3=$compare.gz
    compare=$comp3
fi

tabix -f -p vcf $compare

#comparaison vs témoin
vcf-isec -c $union $compare > $minus

##############
# ANNOTATION #
##############

annot=$RESULTS"/"$prefix"_pindel_soap_prism_minus_"${compname%.*}"_annotated.txt"
echo "snpEff -c $4 -o txt -no-downstream -no-upstream spombe $minus > $annot"
snpEff -c $4 -o txt -no-downstream -no-upstream spombe $minus > $annot


