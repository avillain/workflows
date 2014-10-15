# Adrien Villain (avillain@pasteur.fr)
# Stéphane Descorps-Declère (sdescorp@pasteur.fr)
# 2014 - Variant detection pipeline with a Makefile

#
# -----------------------------------------------------------
# Variables
#

.SECONDARY:
.SILENT:

REF_SEQ?=NC_007793.fa
RESULTS?=./results
R_GROUP?="@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\PU:unit1"
PREF=$(subst .fa,,$(REF_SEQ))
FQ_FILES=$(wildcard *.fastq)
MPIL_FILES=${patsubst %.fastq,%.mpileup,$(FQ_FILES)}
REAL_FILES=${patsubst %.fastq,%.realigned.bam,$(FQ_FILES)}

#
# -----------------------------------------------------------
# Rules
#
.PHONY: all clean index

# Main TARGET
DONADA = @echo "[$@] complete."
NO_OUTPUT= "2>/dev/null || true"
# A bit of color
COL_ON = \033[0;32m
COL_OFF = \033[m

all: $(MPIL_FILES) $(REAL_FILES)

# Preprocessing
$(RESULTS):
	mkdir -p $@

%.fa.fai: %.fa
	@samtools faidx $<

$(PREF).bwt $(PREF).ann $(PREF).amb $(PREF).pac $(PREF).sa : $(REF_SEQ)
	@bwa index -a is -p $(PREF) $< 

$(PREF).dict: $(REF_SEQ)
	@CreateSequenceDictionary REFERENCE=$(REF_SEQ) OUTPUT=$@

# Reads mapping

%.sam: $(PREF).bwt $(FQ_FILES)
	@echo -e "[BWA] Mapping ${COL_ON} $< ${COL_OFF}"
	@bwa mem -M -R $(R_GROUP) $(PREF) $(FQ_FILES) > $@ 

%.fixed.sam: %.sam
	@FixMateInformation INPUT=$< OUTPUT=$@

%.unsorted.bam: %.fixed.sam
	@samtools view -bS -q 2 $< > $@

%.sorted.bam: %.unsorted.bam
	@samtools sort $< $*.sorted

%.stats: %.sorted.bam
	@samtools flagstat $< > $@

%.dedup.bam: %.sorted.bam
	@MarkDuplicates INPUT=$< OUTPUT=$@ M=$*.metrics

%.dedup.bai: %.dedup.bam
	@samtools index $<

%.mpileup: %.dedup.bam $(REF_SEQ) %.dedup.bai %.stats
	@samtools mpileup -f $(REF_SEQ) $< > $@


# GATK
%.list: %.dedup.bam $(REF_SEQ) $(PREF).dict
	@GenomeAnalysisTK -T RealignerTargetCreator -R $(REF_SEQ) -I $< -o $@


%.realigned.bam: %.dedup.bam %.list $(REF_SEQ)
	@GenomeAnalysisTK -T IndelRealigner -R $(REF_SEQ) -I $< -targetIntervals $(word 2,$^) -o $@

#GenomeAnalysisTK -T UnifiedGenotyper -R $refgenome -I $realignedbam -ploidy 1 -glm BOTH -stand_call_conf 10 -stand_emit_conf 5 -o $rawvariants

#GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType SNP -o $rawsnps
#filtexpr="DP<$res || FS > 60.0 || MQ < 40.0"
#GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawsnps --filterExpression "$filtexpr" --filterName $filtname -o $filteredsnps

#GenomeAnalysisTK -T SelectVariants -R $refgenome --variant $filteredsnps -select "vc.isNotFiltered()" --out $filteredcovsnps

#grep '#' "$filteredcovsnps" > "$filteredhetsnps" # header
#grep -v '#' "$filteredcovsnps" | awk '{split($10,geno,":");split(geno[2],ad,","); if(ad[1]*100/(ad[1]+ad[2])<5) print $0}' >> "$filteredhetsnps" # filter on sample/genotype field : 95% of reads must give the same information (e.g. 134 A and 67 T isn't a SNP)

#filtindel="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
#GenomeAnalysisTK -T SelectVariants -R $refgenome -V $rawvariants -selectType INDEL -o $rawindels
#GenomeAnalysisTK -T VariantFiltration -R $refgenome -V $rawindels --filterExpression "$filtindel" --filterName "base" -o $filteredindels





# Others

clean:
	@rm -fr *.mpileup *.list *.dict *.metrics *.stats *.bai *.bam *.sam *.fai *.pac *.ann *.amb *.sa *.bwt

print-%:
	@echo "[ $* ] : $($*)"


