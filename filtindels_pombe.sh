#!/bin/bash

#$ -S /bin/bash

#$ -cwd

#$ -q pf4 

#$ -N pindel_test

#$ -l mem_total=25G

module load snpEff/3.5

#### $1 : prefix $2 : snpeff.config $3 : coverage

#filenames
TMP="./tmp"
RESULTS="./results"
prefix=$1
coverage=$3

### Pindel

pindelsubtracted=$RESULTS"/"$prefix"_pindel_subtracted.vcf"
pindvcf_filt=$RESULTS"/"$prefix"_pindel_filtered.vcf"

pindvcf0=$RESULTS"/"$prefix"_pindel_filt0.vcf"
pindvcf25=$RESULTS"/"$prefix"_pindel_filt25.vcf"
pindvcf50=$RESULTS"/"$prefix"_pindel_filt50.vcf"
pindvcf100=$RESULTS"/"$prefix"_pindel_filt100.vcf"

#filtering
grep "#" $pindelsubtracted > $pindvcf_filt
awk -v cov=$coverage '/^I.*/ {split($10,tab,":"); split(tab[2],ad,","); if(ad[2]>=cov && match(tab[1],/[01]\/1/) && $7=="PASS" && (($1=="III" && ($2>23139 || $2<2440994)) || ($1=="I" && ($2>7618 || $2<5569804)) || ($1=="II" && $2<4532901))) print $0}' $pindelsubtracted >> $pindvcf_filt

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

soapvcf_filt=$RESULTS"/"$prefix"_soap_filtered.vcf"
soapsubtracted=$RESULTS"/"$prefix"_soap_subtracted.vcf"

#filtering
grep "##" $soapsubtracted > $soapvcf_filt
grep "#CHR" $soapsubtracted >> $soapvcf_filt
grep -P "^[I].*PASS.*1/1" $soapsubtracted | awk -v cov=$coverage '{split($8,tab,";"); split(tab[1],tab2,"="); if (tab2[2]>=cov && (($1=="III" && ($2>23139 || $2<2440994)) || ($1=="I" && ($2>7618 || $2<5569804)) || ($1=="II" && $2<4532901))) print $0}' >> $soapvcf_filt

### PRISM (doit être lancé pour chaque chromosome)
prismvcf_filt=$RESULTS"/"$prefix"_prism_filtered.vcf"
prismsubtracted=$RESULTS"/"$prefix"_prism_subtracted.vcf"

awk -v cov=$coverage '/^I.*/ {split($8,geno,";");split(geno[1],ad,"="); if (ad[2]<cov || ($1=="III" && ($2<=23139 || $2>=2440994)) || ($1=="I" && ($2<=7618 || $2>=5569804)) || ($1=="II" && $2>=4532901)); else print $0}' $prismsubtracted > $prismvcf_filt

echo -e "##fileformat=VCFv4.0\n####fileDate=20140724\n####source=prism\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total number of reads in haplotype window">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n###ALT=<ID=ALTER,Description="Alter">\n###FILTER=<ID=q10,Description="Quality below 10">\n###FILTER=<ID=hp10,Description="Reference homopolymer length was longer than 10">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | cat - $prismvcf_filt > dum && mv dum $prismvcf_filt

###############
# COMPARAISON #
###############

fusion_filt=$RESULTS"/"$prefix"_pindel_soap_prism_fusion_filtered_minusPB1623.vcf"
minus_filt=$RESULTS"/"$prefix"_pindel_soap_prism_filtered_minusPB1623.vcf"
mpil="/pasteur/projets/NGS-Dyngen/snps/PB1623/G21-1623-1/tmp/G21-1623-1.mpileup"

#joinx vcf-merge -e -s $pindelsubtracted $soapsubtracted $prismsubtracted > $minus
combinevcf.py $pindvcf_filt $soapvcf_filt $prismvcf_filt $mpil $minus_filt > $fusion_filt

##############
# ANNOTATION #
##############

annot_filt=$RESULTS"/"$prefix"_pindel_soap_prism_filtered_minusPB1623_annotated.txt"
snpEff -c $2 -o txt -no-downstream -no-upstream spombe $minus_filt > $annot_filt


