#!/bin/bash
## avillain@pasteur.fr
## 25/10/13
## To be used on the cluster

## Options to add :
####### - set k value for assembly (+ automatic k determination from reads length)
####### - set taxon id and config file for snp annotation using SnpEff
####### - ...

## Other improvments :
####### - better files/directories structure (log files, etc)
####### - better error handling + restart after failure
####### - write a python script to extract SNP flanking regions

usage()
{
cat << EOF
usage: $0 options
This script will run assembly for 1 strain and use it as a reference for SNP calling of another strain. Only single-end files are supported at the moment.

OPTIONS:
   -h      Show this message
   -q      Query reads (mandatory, .fastq format)
   -r      Reference reads (mandatory, .fastq format)
   -Q      Name for query (default = basename of q)
   -R      Name for reference (default = basename of r)
   -o      Output directory (default = .)
   -k      kmer size for assembly (default = estimated with kmergenie - http://kmergenie.bx.psu.edu/ Chikhi R., Medvedev P. Informed and Automated k-Mer Size Selection for Genome Assembly, HiTSeq 2013)
   -n      Don't perform quality control
   -b      Check reference reads vs reference assembly
   -v      Verbose

EXAMPLE:
./script.sh -q query.fastq -r ref.fastq -Q strainA -R strainB -o results -k 51 -n
EOF
}

QUERY=
REFERENCE=
QNAME=
REFNAME=
OUTPUT=
KMER=
NOQC=
VERBOSE=
while getopts “hq:r:o:Q:R:o:k:nbv” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         q)
             QUERY=$OPTARG
             ;;
         r)
             REFERENCE=$OPTARG
             ;;
         Q)
	     QNAME=$OPTARG
	     ;;
	 R)
             REFNAME=$OPTARG
             ;;
         o)
	     OUTPUT=$OPTARG
	     ;;
	 k)
	     KMER=$OPTARG
	     ;;
	 n)
	     NOQC=1
	     ;;
	 b)
	     BACKWARD=1
	     ;;
	 v)
             VERBOSE=1
             ;;
     ?)
         usage
         exit
         ;;
     esac
done

if [[ -z $QUERY ]] || [[ -z $REFERENCE ]]
then
     echo "You did not supply reference reads (-r) or query reads (-q)"
     usage
     exit 1
fi

#Versions
module load bwa/0.7.5a
module load GenomeAnalysisTK/2.7-2
module load kmergenie

which bwa

#prefixes and directories of outputs
if [[ -z $QNAME ]]
then
     QNAME=$(basename $QUERY)
     QNAME="${QNAME%.*}"
fi

if [[ -z $REFNAME ]]
then
     REFNAME=$(basename $REFERENCE})
     REFNAME="${REFNAME%.*}"

fi

if [[ -z $OUTPUT ]]
then
     OUTPUT='.'
else
     OUTPUT=${OUTPUT%/}
     if [[ ! -s $OUTPUT ]]
     then
          mkdir $OUTPUT
     fi
fi

TMP=$OUTPUT"/tmp"
RESULTS=$OUTPUT"/results"

if [[ ! -s $TMP ]]
then
     mkdir $TMP
fi

if [[ ! -s $RESULTS ]]
then
     mkdir $RESULTS
fi

PREFIX=$QNAME"vs"$REFNAME
readgroup="@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\PU:unit1"

###################
# Quality control #
###################

# fqCleaner must be in the $PATH

if [[ -z $NOQC ]]
then
     REFQC=$RESULTS"/"$(basename ${REFERENCE%.*}"_clean.fastq")
     QRYQC=$RESULTS"/"$(basename ${QUERY%.*}"_clean.fastq")
     fqCleaner.sh -f $REFERENCE -x $REFQC -s "QFAD"
     fqCleaner.sh -f $QUERY -x $QRYQC -s "QFAD"
     REFERENCE=$REFQC
     QUERY=$QRYQC
fi

####################
# De novo assembly #
####################
velvetdir=$TMP"/"$REFNAME"_assembly"

# if assembly has not been done yet
if [[ ! -s $velvetdir"/contigs.fa" ]]
then
     # if k wasn't defined by user, estimate it using kmergenie (http://kmergenie.bx.psu.edu/ Chikhi R., Medvedev P. Informed and Automated k-Mer Size Selection for Genome Assembly, HiTSeq 2013)
     if [[  -z $KMER ]]
     then
          kmerout=$TMP"/out.txt"
          kmergenie -o $TMP"/kmerg/histo" $REFERENCE > $kmerout
          KMER="$(grep 'best k:' $kmerout)"; KMER="${KMER##*:}";
     fi
     velveth-101 $velvetdir $KMER -short -fastq $REFERENCE
     velvetg-101 $velvetdir -cov_cutoff 4 -min_contig_lgth 200
fi

refgenome=$RESULTS"/"$REFNAME".fa"
cp $velvetdir"/contigs.fa" $refgenome

######################
# Strains similarity #
######################

# Optional - not done yet

#############
# Alignment #
#############

# dispatch files in /tmp and /results

outputsam=$TMP"/"$PREFIX".sam"
outputbam=$TMP"/"$PREFIX".bam"
dedupbam=$RESULTS"/"$PREFIX".bam"
dictname=$RESULTS"/"$REFNAME".dict"
metrics=$TMP"/"$REFNAME".metrics"
rawvariants=$RESULTS"/"$PREFIX"_raw.vcf"

logfile=$RESULTS"/"$PREFIX".log"

#Index reference genome
fai=$refgenome".fai"

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

#Mapping

bwa mem -M -R $readgroup $RESULTS"/"$REFNAME $QUERY > $outputsam

#Sam to Bam
SortSam INPUT=$outputsam OUTPUT=$outputbam SORT_ORDER="coordinate"

##Mark Duplicates
MarkDuplicates INPUT=$outputbam OUTPUT=$dedupbam M=$metrics

BuildBamIndex INPUT=$dedupbam

###############
# SNP calling #
###############

GenomeAnalysisTK -T UnifiedGenotyper -R $refgenome -I $dedupbam -ploidy 1 -glm BOTH -stand_call_conf 30 -stand_emit_conf 10 -o $rawvariants

###SNPs
##filenames and filters
rawsnps=$TMP"/"$PREFIX"_snps_raw.vcf"
filtexpr="DP<50 || FS > 60.0 || MQ < 40.0"
filtname="cov50"
filteredsnps=$TMP"/"$PREFIX"snps_cov50.vcf"
filteredcovsnps=$TMP"/"$PREFIX"snps_cov50selected.vcf"
filteredhetsnps=$RESULTS"/"$PREFIX"snps_filtered.vcf"

GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType SNP -o $rawsnps
GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawsnps --filterExpression "$filtexpr" --filterName $filtname -o $filteredsnps
GenomeAnalysisTK -T SelectVariants -R $refgenome --variant $filteredsnps -select "vc.isNotFiltered()" --out $filteredcovsnps

grep '#' "$filteredcovsnps" > "$filteredhetsnps" # header
grep -v '#' "$filteredcovsnps" | awk '{split($10,geno,":");split(geno[2],ad,","); if(ad[1]*100/(ad[1]+ad[2])<5) print $0}' >> "$filteredhetsnps" # filter on sample/genotype field : 95% of reads must give the same information (e.g. 134 A and 67 T isn't a SNP)

#############################
#Mapping of Reference strain#
#############################

if [[ ! -z $BACKWARD ]]
then
	PREFIXBIS=$REFNAME"vs"$REFNAME
	
	outputsam=$TMP"/"$PREFIXBIS".sam"
	outputbam=$TMP"/"$PREFIXBIS".bam"
	dedupbam=$RESULTS"/"$PREFIXBIS".bam"
	dictname=$RESULTS"/"$PREFIXBIS".dict"
	metrics=$TMP"/"$PREFIXBIS".metrics"
	rawvariants=$RESULTS"/"$PREFIXBIS"_raw.vcf"
	
	bwa mem -M -R $readgroup $RESULTS"/"$REFNAME $REFERENCE > $outputsam
	
	#Sam to Bam
	SortSam INPUT=$outputsam OUTPUT=$outputbam SORT_ORDER="coordinate"
	
	##Mark Duplicates
	MarkDuplicates INPUT=$outputbam OUTPUT=$dedupbam M=$metrics
	
	BuildBamIndex INPUT=$dedupbam
	
	###############
	# SNP calling #
	###############
	
	GenomeAnalysisTK -T UnifiedGenotyper -R $refgenome -I $dedupbam -ploidy 1 -glm BOTH -stand_call_conf 30 -stand_emit_conf 10 -o $rawvariants
	
	###SNPs
	##filenames and filters
	rawsnps=$TMP"/"$PREFIXBIS"_snps_raw.vcf"
	filtexpr="DP<50 || FS > 60.0 || MQ < 40.0"
	filtname="cov50"
	filteredsnps=$TMP"/"$PREFIXBIS"snps_cov50.vcf"
	filteredcovsnps=$TMP"/"$PREFIXBIS"snps_cov50selected.vcf"
	filteredhetsnps=$RESULTS"/"$PREFIXBIS"snps_filtered.vcf"
	
	GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType SNP -o $rawsnps
	GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawsnps --filterExpression "$filtexpr" --filterName $filtname -o $filteredsnps
	GenomeAnalysisTK -T SelectVariants -R $refgenome --variant $filteredsnps -select "vc.isNotFiltered()" --out $filteredcovsnps
	
	grep '#' "$filteredcovsnps" > "$filteredhetsnps" # header
	grep -v '#' "$filteredcovsnps" | awk '{split($10,geno,":");split(geno[2],ad,","); if(ad[1]*100/(ad[1]+ad[2])<5) print $0}' >> "$filteredhetsnps" # filter on sample/genotype field : 95% of reads must give the same information (e.g. 134 A and 67 T isn't a SNP)
fi
