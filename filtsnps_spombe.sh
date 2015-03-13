#!/bin/bash

#./filtsnps_spombe.sh J60_01 15 /pasteur/projets/NGS-Dyngen/references/spombe223.fa /pasteur/projets/NGS-Dyngen/references/PB1623/PB1623snps_raw.vcf.gz /pasteur/projets/NGS-Dyngen/references/PB1623/PB1623indels_raw.vcf.gz

PREFIX=$1

TMP="./tmp"
RESULTS="./results"

res=$2

refgenome=$3

if [ ! -f ${refgenome%.*}.fai ]
then
        samtools faidx $refgenome
fi

if [ ! -f ${refgenome%.*}.dict ]
then
        CreateSequenceDictionary REFERENCE=$refgenome OUTPUT=${refgenome%.*}.dict
fi

SNPSUB=$4

INDSUB=$5

echo "[info] Minimum coverage set to $res"

###SNPs
##filenames and filters
rawsnps=$TMP"/"$PREFIX"_snps_raw.vcf"
filtexpr="DP<$res || FS > 60.0 || MQ < 40.0"
filtname="cov"$res
filteredsnps=$TMP"/"$PREFIX"snps_cov"$res".vcf"
#filteredcovsnps=$TMP"/"$PREFIX"snps_cov"$res"selected.vcf"
filteredhetall=$TMP"/"$PREFIX"snpsall_filtered"$res".vcf"
filteredhetsnps=$RESULTS"/"$PREFIX"snps_filtered"$res".vcf"

echo "[info] filtering SNPs using expression $filtexpr"
echo '[cmd] GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType SNP -o $rawsnps'
GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawsnps --filterExpression "$filtexpr" --filterName $filtname -o $filteredsnps 
GenomeAnalysisTK -T SelectVariants -R $refgenome --variant $filteredsnps -select "vc.isNotFiltered()" --out $filteredhetall 

grep "#" $filteredhetall > $filteredhetsnps
awk '/^I.*/ { if (($1=="III" && ($2<=23139 || $2>=2440994)) || ($1=="I" && ($2<=7618 || $2>=5569804)) || ($1=="II" && $2>=4532901)); else print $0}' $filteredhetall >> $filteredhetsnps

###Indels
rawindels=$TMP"/"$PREFIX"_indels_raw.vcf"
filteredindels_tmp=$TMP"/"$PREFIX"indels_filtering"$res".vcf"
filteredallindels=$TMP"/"$PREFIX"indelsall_filtered"$res".vcf"
filteredindels=$RESULTS"/"$PREFIX"indels_filtered"$res".vcf"
filtindel="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawindels --filterExpression "$filtindel" --filterName "base" -o $filteredindels_tmp 
GenomeAnalysisTK -T SelectVariants -R $refgenome --variant $filteredindels_tmp -select "vc.isNotFiltered()" --out $filteredallindels

grep "#" $filteredallindels > $filteredindels
awk '/^I.*/ { if (($1=="III" && ($2<=23139 || $2>=2440994)) || ($1=="I" && ($2<=7618 || $2>=5569804)) || ($1=="II" && $2>=4532901)); else print $0}' $filteredallindels >> $filteredindels

if [[ ! -z $SNPSUB ]] && [[ ! -z $INDSUB ]]
then
	echo "[info] substracting control SNPs and indels from results"
	snpminus=$RESULTS"/"$PREFIX"snps_minus_union1623_"$res".vcf"
	indelminus=$RESULTS"/"$PREFIX"indels_minus_union1623_"$res".vcf"
	bgzip $filteredhetsnps
	tabix -p vcf $filteredhetsnps.gz
	vcf-isec -c $filteredhetsnps.gz $SNPSUB > $snpminus
	gunzip $filteredhetsnps.gz
	rm $filteredhetsnps.gz.tbi
	bgzip $filteredindels
	tabix -p vcf $filteredindels.gz
	vcf-isec -c $filteredindels.gz $INDSUB > $indelminus
	gunzip $filteredindels.gz
	rm $filteredindels.gz.tbi
	allvariants=$RESULTS"/"$PREFIX"allvariants_filtered_minus_union1623_"$res".vcf"
	tmpvar=$TMP"/allvariants_minus1623_tmp_"$res".vcf"
	snpeffannotated=$RESULTS"/"$PREFIX"allvariants_annotated_minus_union1623_"$res".vcf"
else
	echo "[info] no substraction of control SNPs and indels from results"
	snpminus=$filteredhetsnps
	indelminus=$filteredindels
	allvariants=$RESULTS"/"$PREFIX"allvariants_filtered_"$res".vcf"
	tmpvar=$TMP"/allvariants.vcf"
	snpeffannotated=$RESULTS"/"$PREFIX"allvariants_annotated_"$res".vcf"
fi

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
        snpstats=$TMP"/"$PREFIX"annot"
        echo "[cmd] snpEff -c $configfile -v -s $snpstats -o gatk $refID $allvariants > $snpeffannotation"
        snpEff -c $configfile -v -s $snpstats -o gatk $refID $allvariants > $snpeffannotation
        GenomeAnalysisTK -T VariantAnnotator -R $refgenome -A SnpEff --variant $allvariants --snpEffFile $snpeffannotation -o $snpeffannotated
fi

