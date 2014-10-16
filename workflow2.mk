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
FILT_EXPR="DP<50 || FS > 60.0 || MQ < 40.0"
FILT_NAME="cov50"
FILT_INDEL="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
PREF=$(subst .fa,,$(REF_SEQ))
FQ_FILES=$(wildcard *.fastq)
MPIL_FILES=${patsubst %.fastq,%.mpileup,$(FQ_FILES)}
REAL_FILES=${patsubst %.fastq,%.realigned.bam,$(FQ_FILES)}
VCF_FILES=${patsubst %.fastq,%.filtered.vcf,$(FQ_FILES)}

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

all: $(MPIL_FILES) $(VCF_FILES)

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

%.raw.vcf: %.realigned.bam $(REF_SEQ)
	@GenomeAnalysisTK -T UnifiedGenotyper -R $(REF_SEQ) -I $< -ploidy 1 -glm BOTH -stand_call_conf 10 -stand_emit_conf 5 -o $@

%.raw.snp.vcf: %.raw.vcf $(REF_SEQ)
	@GenomeAnalysisTK -T SelectVariants -R $(REF_SEQ) -V $< -selectType SNP -o $@

%.filtering.snp.vcf: %.raw.snp.vcf $(REF_SEQ)
	@GenomeAnalysisTK -T VariantFiltration -R $(REF_SEQ) -V $< --filterExpression $(FILT_EXPR) --filterName $(FILT_NAME) -o $@

%.filtered_1.snp.vcf: %.filtering.snp.vcf $(REF_SEQ)
	@GenomeAnalysisTK -T SelectVariants -R $(REF_SEQ) --variant $< -select "vc.isNotFiltered()" --out $@

%.filtered_2.snp.vcf: %.filtered_1.snp.vcf 
	@grep '#' $< > $@ # header
	@grep -v '#' $< | awk '{split($$10,geno,":");split(geno[2],ad,","); if(ad[1]*100/(ad[1]+ad[2])<5) print $$0}' >> $@

%.raw.indel.vcf: %.raw.vcf $(REF_SEQ)
	@GenomeAnalysisTK -T SelectVariants -R $(REF_SEQ) -V $< -selectType INDEL -o $@

%.filtering.indel.vcf: %.raw.indel.vcf $(REF_SEQ)
	@GenomeAnalysisTK -T VariantFiltration -R $(REF_SEQ) -V $< --filterExpression $(FILT_INDEL) --filterName "base" -o $@

%.filtered.indel.vcf: %.filtering.indel.vcf $(REF_SEQ)
	@GenomeAnalysisTK -T SelectVariants -R $(REF_SEQ) --variant $< -select "vc.isNotFiltered()" --out $@

%.filtered.vcf: %.raw.vcf %.filtered_2.snp.vcf %.filtered.indel.vcf
	@grep '#' $< > $@
	@grep -v '#' $(word 2,$^) $(word 3,$^) |sort -k2n,2 >> $@
	


# Others

clean:
	@rm -fr *.vcf *.vcf.idx *.mpileup *.list *.dict *.metrics *.stats *.bai *.bam *.sam *.fai *.pac *.ann *.amb *.sa *.bwt

print-%:
	@echo "[ $* ] : $($*)"


