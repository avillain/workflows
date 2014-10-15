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

cp $2 $reference
cp $3 $compare

#insert size
insertfile=$TMP"/"$prefix"_insertsize.txt"

if [ ! -f $insertfile ]
then
    $HOME/scripts/getinsertsize.py $1 > $insertfile
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
pindel8cols=$RESULTS"/"$prefix"_pindel_filtered_8cols.vcf"

#configuration
echo -e "$sortedbamfile\t$insertmean\t$prefix" > $pindelconfigfile

#detection
pindel -f $reference -i $pindelconfigfile -c ALL -o $TMP"/"$prefix

#conversion
pindel2vcf -p $pindout -r $reference -R spombe -d 09052011 -v $pindvcf

#filtering
grep "#" $pindvcf > $pindvcf_filt
awk '!/^#.*/ {split($10,tab,","); if(tab[2]>=5) print $0}' $pindvcf >> $pindvcf_filt

grep "^##" $pindvcf_filt > $pindel8cols
awk 'BEGIN {OFS = "\t"} !/^##.*/ { print $1, $2, $3, $4, $5, $6, $7, $8 }' $pindvcf_filt >> $pindel8cols

### SOAPindel

soapvcf=$TMP"/"$prefix"_soap.vcf"
soapvcf_filt=$RESULTS"/"$prefix"_soap_filtered.vcf"
soap8cols=$RESULTS"/"$prefix"_soap_filtered_8cols.vcf"

mappinglist=$TMP"/"$prefix"_mapping.list"
readslist=$TMP"/"$prefix"_reads.list"
searchprefix=$prefix
#searchprefix="1623-1"

#configuration
ls /pasteur/projets/NGS-Dyngen/fastq_files/raw_files/*/*$searchprefix*.fastq > $readslist
echo -e "$sortedbamfile\t$insertmean\t$insertsd\t$readslength\tPAIR" > $mappinglist

#detection
indel_detection.ibam.pl $mappinglist $reference $readslist -wd $TMP

#merging chromosomes files
grep "#" $TMP/result/*chr1*/*indel.vcf > $soapvcf
grep -hv "#" $TMP/result/*/*indel.vcf >> $soapvcf

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
for i in $refdir; do chrname=`basename ${i%.*}`; j=$TMP"/"$prefix"_"$chrname".sam"; grep "SN:$chrname" $1 > $j; head -n 7 $1 | tail -n 1 >> $j; grep -P "\t$chrname\t" $1 >> $j; run_PRISM.sh -m $insertmean -e $insertsd -r $i -i $j -I $TMP -O $TMP"/"$chrname; done

#filtering
for i in $refdir; do chrname=`basename ${i%.*}`; awk '{ if ($6>=5) print $0 }' $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv" > $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv_supported"; done

#conversion
for i in $refdir; do chrname=`basename ${i%.*}`; out=$TMP"/"$prefix"_"$chrname"_prism.vcf"; prism2vcf.py -f $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv_supported" -o $out; done
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
comparegz=$compare.gz
uniongz=$union.gz
compname=`basename $compare`
minus=$RESULTS"/"$prefix"_pindel_soap_prism_minus_"${compname%.*}".vcf"

bgzip $pindel8cols
bgzip $soap8cols
bgzip $prism8cols
tabix -p vcf $pindgz
tabix -p vcf $soapgz
tabix -p vcf $prismgz

#combinaison des outils
vcf-isec -a -n +1 $pindgz $soapgz $prismgz > $union

if [ ! -f $uniongz ]
then
    bgzip $union
fi

tabix -f -p vcf $uniongz

if [ ! -f $comparegz ]
then
    bgzip $compare
fi

tabix -f -p vcf $comparegz

#comparaison vs témoin
vcf-isec -c $uniongz $comparegz > $minus

##############
# ANNOTATION #
##############

annot=$RESULTS"/"$prefix"_pindel_soap_prism_minus_"${compname%.*}"_annotated.txt"
echo "snpEff -c $4 -o txt -no-downstream -no-upstream spombe $minus > $annot"
snpEff -c $4 -o txt -no-downstream -no-upstream spombe $minus > $annot


