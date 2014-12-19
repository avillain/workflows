#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script run the SNP and Indels analysis for one genome using fqCleaner for reads pre-processing (optional), bwa mem for reads mapping, GATK2 for variant calling, and snpEff for variant annotation (optional)

OPTIONS:
   -h      Show this message
   -i      Input reads (.fastq format)
   -j	   Input mate reads (.fastq format)
   -r      Reference genome (.fasta format)
   -p      Output files prefix
   -n      No quality control performed
   -a      Annotation of variants with SnpEff (see -c and -d options)
   -c	   SnpEff config file (mandatory if -a)
   -d      SnpEff database ID of the species  (mandatory if -a)
   -v      Verbose

EXAMPLE:
./script.sh -i reads1.fq -j reads2.fq -r ref.fasta -p prefix -n -a -c SNPeffconfig1.txt -d spombe
EOF
}

INPUT=
INPUTMATE=
REFERENCE=
VERBOSE=
NOQC=
ANNOT=
CONF=
ID=
while getopts “hi:j:r:o:p:c:d:vna” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             INPUT=$OPTARG
             ;;
	 j)
	     INPUTMATE=$OPTARG
	     ;;
         r)
             REFERENCE=$OPTARG
             ;;
         p)
             PREFIX=$OPTARG
             ;;
         c)
             CONF=$OPTARG
             ;;
         d)
             ID=$OPTARG
             ;;
         v)
             VERBOSE=1
             ;;
         n)
             NOQC=1
             ;;
         a)
             ANNOT=1
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

echo "[info] parsing command-line"

if [[ -z $INPUT ]] || [[ -z $REFERENCE ]]
then
     echo "ERROR : Please supply a .fasta reference genome (-r) and .fastq reads to compare (-i)"
     usage
     exit 1
fi

if [[ -z $INPUTMATE ]]
then
	echo "[info] no mate reads file supplied, input is single reads"
fi

if [[ ! -z $ANNOT ]] && ( [[ -z $CONF ]] || [[ -z $ID ]] )
then
	echo "ERROR : SnpEff config file (-c) and SnpEff database ID (-d) are mandatory to perform SnpEff annotation (-a)"
	usage
	exit 1
fi

#basenames
refname=$(basename "$REFERENCE")
extension="${refname##*.}"
refname="${refname%.*}"
readgroup="@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\PU:unit1"

#prefix of outputs
if [[ -z $PREFIX ]]
then
      dir=$(dirname $INPUT)
      PREFIX=${name%.*}"vs"$refname 
      echo "[info] no prefix supplied, output files will begin with $PREFIX"
else
      dir=$(dirname $PREFIX)
      PREFIX=$(basename $PREFIX)
      echo "[info] output files will begin with $PREFIX"
fi

TMP=$dir"/tmp"
RESULTS=$dir"/results"
echo "[info] working in $dir/tmp - results stored in $dir/results"
if [[ ! -s $RESULTS ]]
     then
          mkdir $RESULTS
fi

if [[ ! -s $TMP ]]
     then
          mkdir $TMP
fi

outputsam=$TMP"/"$PREFIX".sam"
fixedsam=$TMP"/"$PREFIX"_fixed.sam"
outputbam=$TMP"/"$PREFIX".bam"
dedupbam=$TMP"/dedup_"$PREFIX".bam"
dictname=$RESULTS"/"$refname".dict"
metrics=$TMP"/"$refname".metrics"
realignedbam=$TMP"/realigned_"$PREFIX".bam"
rawvariants=$TMP"/"$PREFIX"_raw.vcf"

logfile=$TMP"/"$PREFIX".log"
REFNAME=$(basename $REFERENCE})
REFNAME="${REFNAME%.*}"
refgenome=$RESULTS"/"$REFNAME".fa"

if [ ! -f $refgenome ]
then 
      cp $REFERENCE $refgenome
fi

#Versions
echo "[info] versions : bwa/0.7.5a GenomeAnalysisTK/2.7-2 snpEff/3.5"
module load bwa/0.7.5a
module load GenomeAnalysisTK/2.7-2
module load snpEff/3.5

#Index reference genome
fai=$TMP"/"$REFERENCE".fai"
echo "[info] indexing reference genome"
if [ ! -f $fai ]
then
	samtools faidx $refgenome 
fi

if [ ! -f $RESULTS"/"$REFNAME".bwt" ]
then
        bwa index -a is -p $RESULTS"/"$REFNAME $refgenome 
fi

if [ ! -f $dictname ]
then
        CreateSequenceDictionary REFERENCE=$refgenome OUTPUT=$dictname 
fi

###################
# Quality control #
###################
# fqCleaner must be in the $PATH

if [[ -z $NOQC ]]
then
     INQC=$RESULTS"/"$(basename ${INPUT%.*}"_clean.fastq")
     if [[ -z $INPUTMATE ]]
     then
          echo "[info] performing single-end quality control using fqCleaner"
          echo '[cmd] fqCleaner.sh -f $INPUT -x $INQC -s "QFAD"'
          fqCleaner.sh -f $INPUT -x $INQC -s "QFAD" -l 80 -q 30 || exit 1
     else
	  INQCMATE=$RESULTS"/"$(basename ${INPUTMATE%.*})"_clean.fastq"
          echo "[info] performing paired-end quality control using fqCleaner"
          echo "[cmd] fqCleaner.sh -f $INPUT -r $INPUTMATE -x $INQC -y $INQCMATE -s "QFAD" -l 80 -q 30 "
          fqCleaner.sh -f $INPUT -r $INPUTMATE -x $INQC -y $INQCMATE -s "QFAD" -l 80 -q 30 || exit 1
     fi
     INPUT=$INQC
     INPUTMATE=$INQCMATE
fi

if [[ ! -z $INPUTMATE ]]
then
     echo '[cmd] bwa mem -M -R '$readgroup $RESULTS"/"$REFNAME $INPUT $INPUTMATE' > $outputsam'
     bwa mem -M -R $readgroup $RESULTS"/"$REFNAME $INPUT $INPUTMATE > $outputsam || exit 1
else
     echo '[cmd] bwa mem -M -R '$readgroup $RESULTS"/"$REFNAME $INPUT' > $outputsam'
     bwa mem -M -R $readgroup $RESULTS"/"$REFNAME $INPUT > $outputsam || exit 1
fi

#Sam to Bam
echo '[cmd] FixMateInformation INPUT=$outputsam OUTPUT=$fixedsam'
FixMateInformation INPUT=$outputsam OUTPUT=$fixedsam

#echo '[cmd] SortSam INPUT=$outputsam OUTPUT=$outputbam SORT_ORDER="coordinate"'
SortSam INPUT=$fixedsam OUTPUT=$outputbam SORT_ORDER="coordinate" 
echo '[info] get mapping statistics'
echo "[cmd] samtools flagstat $outputbam"
samtools flagstat $outputbam


##Mark Duplicates
echo '[cmd] MarkDuplicates INPUT=$outputbam OUTPUT=$dedupbam M=$metrics'
MarkDuplicates INPUT=$outputbam OUTPUT=$dedupbam M=$metrics 

echo '[cmd] BuildBamIndex INPUT=$dedupbam'
BuildBamIndex INPUT=$dedupbam 

#Indel realignment
echo "[info] indels realignment"
echo '[cmd] GenomeAnalysisTK -T RealignerTargetCreator -R $refgenome -I $dedupbam -o $TMP"/"$PREFIX"_target_intervals.list"'
GenomeAnalysisTK -T RealignerTargetCreator -R $refgenome -I $dedupbam -o $TMP"/"$PREFIX"_target_intervals.list" 

echo '[cmd] GenomeAnalysisTK -T IndelRealigner -R $refgenome -I $dedupbam -targetIntervals $TMP"/"$PREFIX"_target_intervals.list" -o $realignedbam'
GenomeAnalysisTK -T IndelRealigner -R $refgenome -I $dedupbam -targetIntervals $TMP"/"$PREFIX"_target_intervals.list" -o $realignedbam 

pileup=$TMP"/"$PREFIX".mpileup"
cover=$TMP"/"$PREFIX".txt"
covplot=$TMP"/"$PREFIX"_coverage.jpeg"

#adjust coverage limit for snp filtering
echo "[info] adjusting coverage limit for snp filtering"
echo '[cmd] samtools mpileup -f $refgenome $realignedbam > $pileup'
#echo '[cmd] cov=`awk '{ total += $4;count++ } END {print total/count}' $pileup`'
#echo '[cmd] res=$(echo "$cov / 2" |bc )'
samtools mpileup -f $refgenome $realignedbam > $pileup 
#cov=`awk '{ total += $4;count++ } END {print total/count}' $pileup`
#res=$(echo "$cov / 2" |bc )
res=15
#coverage plot
#echo "[info] plot mapping coverage along genome"
#echo "[cmd] cat $pileup | cut -f 2,4 > $cover | awk + gnuplot magic"

#cat $pileup | cut -f 2,4 > $cover

#awk 'BEGIN { OFS = "\t" } ;($1-p2)>1{
#for(i=p2+1;i<$1;i++)
#print i,0
#}
#{p2=$1}1' $cover | gnuplot -e "set term jpeg ; set output '$covplot' ; set ylabel 'Position on genome' ; set xlabel 'Coverage' ;  set style data linespoints ; set title 'Coverage along reference' ; set key reverse Left outside ; set grid ; plot '<cat' using 1:2 lt 1 lw 2 smooth bezier title '$PREFIX'"

#Calling
##GATK
echo '[info] SNPs and Indels calling'
echo '[cmd] GenomeAnalysisTK -T UnifiedGenotyper -R $refgenome -I $realignedbam -ploidy 1 -glm BOTH -stand_call_conf 10 -stand_emit_conf 5 -o $rawvariants'
GenomeAnalysisTK -T UnifiedGenotyper -R $refgenome -I $realignedbam -ploidy 1 -glm BOTH -stand_call_conf 10 -stand_emit_conf 5 -o $rawvariants 

###SNPs
##filenames and filters
rawsnps=$TMP"/"$PREFIX"_snps_raw.vcf"
filtexpr="DP<$res || FS > 60.0 || MQ < 40.0"
filtname="cov"$res
filteredsnps=$TMP"/"$PREFIX"snps_cov"$res".vcf"
#filteredcovsnps=$TMP"/"$PREFIX"snps_cov"$res"selected.vcf"
filteredhetall=$TMP"/"$PREFIX"snpsall_filtered.vcf"
filteredhetsnps=$RESULTS"/"$PREFIX"snps_filtered.vcf"

echo "[info] filtering SNPs using expression $filtexpr"
echo '[cmd] GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType SNP -o $rawsnps'

GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType SNP -o $rawsnps 
GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawsnps --filterExpression "$filtexpr" --filterName $filtname -o $filteredsnps 
GenomeAnalysisTK -T SelectVariants -R $refgenome --variant $filteredsnps -select "vc.isNotFiltered()" --out $filteredhetall 

grep "#" $filteredhetall > $filteredhetsnps
awk '/^I.*/ { if (($1=="III" && ($2<=23139 || $2>=2440994)) || ($1=="I" && ($2<=7618 || $2>=5569804)) || ($1=="II" && $2>=4532901)); else print $0}' $filteredhetall >> $filteredhetsnps

#echo '[cmd] grep "#" "$filteredhetsnps" > "$filteredhetsnps"'
#echo "[cmd] grep -v \"#\" $filteredcovsnps | awk '{split($10,geno,":");split(geno[2],ad,\",\"); if(ad[1]*100/(ad[1]+ad[2])<5) print $0}' >> \"$filteredhetsnps\""
#grep '#' "$filteredcovsnps" > "$filteredhetsnps" # header
#grep -v '#' "$filteredcovsnps" | awk '{split($10,geno,":");split(geno[2],ad,","); if(ad[1]*100/(ad[1]+ad[2])<5) print $0}' >> "$filteredhetsnps" # filter on sample/genotype field : 95% of reads must give the same information (e.g. 134 A and 67 T isn't a SNP)

###Indels
rawindels=$TMP"/"$PREFIX"_indels_raw.vcf"
filteredindels_tmp=$TMP"/"$PREFIX"indels_filtering.vcf"
filteredallindels=$TMP"/"$PREFIX"indelsall_filtered.vcf"
filteredindels=$RESULTS"/"$PREFIX"indels_filtered.vcf"
filtindel="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType INDEL -o $rawindels 
GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawindels --filterExpression "$filtindel" --filterName "base" -o $filteredindels_tmp 
GenomeAnalysisTK -T SelectVariants -R $refgenome --variant $filteredindels_tmp -select "vc.isNotFiltered()" --out $filteredallindels

grep "#" $filteredallindels > $filteredindels
awk '/^I.*/ { if (($1=="III" && ($2<=23139 || $2>=2440994)) || ($1=="I" && ($2<=7618 || $2>=5569804)) || ($1=="II" && $2>=4532901)); else print $0}' $filteredallindels >> $filteredindels

snpminus=$RESULTS"/"$PREFIX"snps_minus_union1623.vcf"
indelminus=$RESULTS"/"$PREFIX"indels_minus_union1623.vcf"

bgzip $filteredhetsnps
tabix -p vcf $filteredhetsnps.gz
vcf-isec -c $filteredhetsnps.gz /pasteur/projets/NGS-Dyngen/snps/PB1623/union/PB1623snps_raw.vcf.gz > $snpminus
gunzip $filteredhetsnps.gz
rm $filteredhetsnps.gz.tbi

bgzip $filteredindels
tabix -p vcf $filteredindels.gz
vcf-isec -c $filteredindels.gz /pasteur/projets/NGS-Dyngen/snps/PB1623/union/PB1623indels_raw.vcf.gz > $indelminus
gunzip $filteredindels.gz
rm $filteredindels.gz.tbi

allvariants=$RESULTS"/"$PREFIX"allvariants_filtered_minus_union1623.vcf"
tmpvar=$TMP"/allvariants_minus1623_tmp.vcf"

echo "[info] merging snps and indels in the same .vcf file"
grep -v '#' $snpminus > $tmpvar
grep -v '#' $indelminus >> $tmpvar

grep '#' $snpminus > $allvariants
grep -P "^MT\t" $tmpvar | sort -k2,2n >> $allvariants
grep -P "^AB325691\t" $tmpvar | sort -k2,2n >> $allvariants
grep -P "^MTR\t" $tmpvar | sort -k2,2n >> $allvariants
grep -P "^III\t" $tmpvar | sort -k2,2n >> $allvariants
grep -P "^II\t" $tmpvar | sort -k2,2n >> $allvariants
grep -P "^I\t" $tmpvar | sort -k2,2n >> $allvariants

#Annotation using snpEff
if [[ ! -z $ANNOT ]]
then
        configfile=$CONF
        refID=$ID
        echo "[info] annotation using snpEff : config file $CONF and ref ID $ID"
        #snps
        snpeffannotation=$TMP"/"$PREFIX"_effects.vcf"
        snpeffannotated=$RESULTS"/"$PREFIX"allvariants_annotated_minus_union1623.vcf"
        snpstats=$TMP"/"$PREFIX"annot"
        echo "[cmd] snpEff -c $configfile -v -s $snpstats -o gatk $refID $allvariants > $snpeffannotation"
        snpEff -c $configfile -v -s $snpstats -o gatk $refID $allvariants > $snpeffannotation
        GenomeAnalysisTK -T VariantAnnotator -R $refgenome -A SnpEff --variant $allvariants --snpEffFile $snpeffannotation -o $snpeffannotated
fi

#if [[ ! -z $ANNOT ]]
#then
#	configfile=$CONF
#	refID=$ID
#	echo "[info] annotation using snpEff : config file $CONF and ref ID $ID"
#	#snps
#	snpeffannotation=$TMP"/"$PREFIX"snps_effects.vcf"
#	snpeffannotated=${filteredhetsnps%.*}"_annotated.vcf"
#	snpstats=$TMP"/"$PREFIX"snpannot"
#	echo "[cmd] snpEff -c $configfile -v -s $snpstats -o gatk $refID $filteredhetsnps > $snpeffannotation"
#	snpEff -c $configfile -v -s $snpstats -o gatk $refID $filteredhetsnps > $snpeffannotation
#	GenomeAnalysisTK -T VariantAnnotator -R $refgenome -A SnpEff --variant $filteredhetsnps --snpEffFile $snpeffannotation  -o $snpeffannotated
#	#indels
#	indeleffannotation=$TMP"/"$PREFIX"indels_effects.vcf"
#	indeleffannotated=${filteredindels%.*}"_annotated.vcf"
#	indelstats=$TMP"/"$PREFIX"indelsannot"
#	echo "[cmd] snpEff -c $configfile -v -s $indelstats -o gatk $refID $filteredindels > $indeleffannotation"
#	snpEff -c $configfile -v -s $indelstats -o gatk $refID $filteredindels > $indeleffannotation
#	GenomeAnalysisTK -T VariantAnnotator -R $refgenome -A SnpEff --variant $filteredindels --snpEffFile $indeleffannotation  -o $indeleffannotated
#	
#	allvariants=$RESULTS"/"$PREFIX"allvariants_annotatedi_minus_union1623.vcf"
#	echo "[info] merging snps and indels in the same .vcf file"
#	cat $snpeffannotated > $allvariants
#	grep -v '#' $indeleffannotated >> $allvariants
#fi

