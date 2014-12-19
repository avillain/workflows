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
module load snpEff/3.5

####################
# PRECONFIGURATION #
####################

#qsub script.sh file.sam reference.fa unionPB1623.vcf snpEff.config PB1623_pindel.vcf.gz PB1623_soap.vcf.gz PB1623_prism.vcf.gz

#### $1 : samfile.sam $2 : reference.fa $3 : compare.vcf.gz $4 : snpeff.config $5 : comparepindel.vcf.gz $6 : comparesoap.vcf.gz compareprism.vcf.gz

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
compare=$TMP"/"`basename $3`

compindel=$5
cp $2 $reference
ln $3 $compare

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
pindvcf_filtgz=$pindvcf_filt.gz
pindelsubtracted=$RESULTS"/"$prefix"_pindel_subtracted.vcf"
pindvcf0=$RESULTS"/"$prefix"_pindel_filt0.vcf"
pindvcf25=$RESULTS"/"$prefix"_pindel_filt25.vcf"
pindvcf50=$RESULTS"/"$prefix"_pindel_filt50.vcf"
pindvcf100=$RESULTS"/"$prefix"_pindel_filt100.vcf"

#configuration
echo -e "$sortedbamfile\t$insertmean\t$prefix" > $pindelconfigfile

#detection
pindel -f $reference -i $pindelconfigfile -c ALL -o $TMP"/"$prefix

#conversion
pindel2vcf -p $pindout -r $reference -R spombe -d 09052011 -v $pindvcf

#filtering
grep "#" $pindvcf > $pindvcf_filt
awk '/^I.*/ {split($10,tab,":"); split(tab[2],ad,","); if(ad[2]>=30 && match(tab[1],/[01]\/1/) && $7=="PASS" && (($1=="III" && ($2>23139 || $2<2440994)) || ($1=="I" && ($2>7618 || $2<5569804)) || ($1=="II" && $2<4532901))) print $0}' $pindvcf >> $pindvcf_filt

#subtraction

bgzip $pindvcf_filt
tabix -p vcf $pindvcf_filtgz

vcf-isec -f -c $pindvcf_filtgz $5 > $pindelsubtracted

gunzip $pindvcf_filtgz
rm $pindvcf_filtgz.tbi

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
soapvcf_filtgz=$soapvcf_filt.gz
soapsubtracted=$RESULTS"/"$prefix"_soap_subtracted.vcf"
mappinglist=$TMP"/"$prefix"_mapping.list"
readslist=$TMP"/"$prefix"_reads.list"
searchprefix=$prefix

#configuration
ls /pasteur/projets/NGS-Dyngen/fastq_files/clean_files/*$searchprefix*.fastq > $readslist
echo -e "$sortedbamfile\t$insertmean\t$insertsd\t$readslength\tPAIR" > $mappinglist

#detection
indel_detection.ibam.pl $mappinglist $reference $readslist -wd $TMP

#merging chromosomes files
grep "#" $TMP/result/I\:*/*indel.vcf > $soapvcf
grep -hv "#" $TMP/result/*/*indel.vcf >> $soapvcf

#filtering
grep "##" $soapvcf > $soapvcf_filt
grep "#CHR" $soapvcf >> $soapvcf_filt
grep -P "^[I].*PASS.*1/1" $soapvcf | awk '{split($8,tab,";"); split(tab[1],tab2,"="); if (tab2[2]>=30 && (($1=="III" && ($2>23139 || $2<2440994)) || ($1=="I" && ($2>7618 || $2<5569804)) || ($1=="II" && $2<4532901))) print $0}' >> $soapvcf_filt

bgzip $soapvcf_filt
tabix -p vcf $soapvcf_filtgz

vcf-isec -f -c $soapvcf_filtgz $6 > $soapsubtracted

gunzip $soapvcf_filtgz
rm $soapvcf_filtgz.tbi

### PRISM (doit être lancé pour chaque chromosome)
prismvcf=$TMP"/"$prefix"_prism_unfilt.vcf"
prismvcf_filt=$RESULTS"/"$prefix"_prism_filtered.vcf"
prismvcf_filtgz=$prismvcf_filt.gz
prismsubtracted=$RESULTS"/"$prefix"_prism_subtracted.vcf"
refdir=$reference".DB.SPLIT/*"

#detection
for i in $refdir; do chrname=`basename ${i%.*}`; run_PRISM.sh -m $insertmean -e $insertsd -r $i -i $1 -I $TMP -O $TMP"/"$chrname; done

#filtering
for i in $refdir; do chrname=`basename ${i%.*}`; awk '{ if ($6>=5) print $0 }' $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv" > $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv_supported"; done

#conversion
for i in $refdir; do chrname=`basename ${i%.*}`; out=$TMP"/"$prefix"_"$chrname"_prism.vcf"; prism2vcf.py $TMP"/"$chrname"/split_all.sam_ns_rmmul_cigar_sorted_sv_supported" $reference $out; done
cat $TMP"/"*prism.vcf > $prismvcf

echo -e "##fileformat=VCFv4.0\n####fileDate=20140724\n####source=prism\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total number of reads in haplotype window">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n###ALT=<ID=ALTER,Description="Alter">\n###FILTER=<ID=q10,Description="Quality below 10">\n###FILTER=<ID=hp10,Description="Reference homopolymer length was longer than 10">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > $prismvcf_filt


#echo -e "##fileformat=VCFv4.0\n###fileDate=20140724\n###source=prism\n###reference=/pasteur/projets/NGS-Dyngen/large_indels/test/J64-G15/./tmp/pombe_09052011.fasta.DB.SPLIT\n###INFO=<ID=PROGRAM,Number=1,Type=String,Description="Total number of reads in haplotype window">\n###INFO=<ID=SVTYPE,Number=1,Type=String,Description="Total number of reads 2 in haplotype window">\n###INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Left flank length and right flank length">\n###FILTER=<ID=q10,Description="Quality below 10">\n###FILTER=<ID=hp10,Description="Reference homopolymer length was longer than 10">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > $prismvcf_filt

awk '/^I.*/ {split($8,geno,";");split(geno[1],ad,"="); if (ad[2]<25 || ($1=="III" && ($2<=23139 || $2>=2440994)) || ($1=="I" && ($2<=7618 || $2>=5569804)) || ($1=="II" && $2>=4532901)); else print $0}' $prismvcf >> $prismvcf_filt

bgzip $prismvcf_filt
tabix -p vcf $prismvcf_filtgz

vcf-isec -f -c $prismvcf_filtgz $7 > $prismsubtracted

gunzip $prismvcf_filtgz
rm $prismvcf_filtgz.tbi


###############
# COMPARAISON #
###############

minus=$RESULTS"/"$prefix"_pindel_soap_prism_minusPB1623.vcf"

joinx vcf-merge -e -s $pindelsubtracted $soapsubtracted $prismsubtracted > $minus

##############
# ANNOTATION #
##############

annot=$RESULTS"/"$prefix"_pindel_soap_prism_minusPB1623_annotated.txt"
echo "snpEff -c $4 -o txt -no-downstream -no-upstream spombe $minus > $annot"
snpEff -c $4 -o txt -no-downstream -no-upstream spombe $minus > $annot


