00.K-mer.sh
## soapec.sh The genome size was predicted based on different K-mer and second-sequencing data
``` bash
KmerFreq_HA -k 21 -t 50 -p PigNoseTurtle -l zb.R1.fq.gz zb.R2.fq.gz -L 150 
KmerFreq_HA -k 23 -t 50 -p PigNoseTurtle -l zb.R1.fq.gz zb.R2.fq.gz -L 150 
KmerFreq_HA -k 25 -t 50 -p PigNoseTurtle -l zb.R1.fq.gz zb.R2.fq.gz -L 150 
KmerFreq_HA -k 27 -t 50 -p PigNoseTurtle -l zb.R1.fq.gz zb.R2.fq.gz -L 150
```
#01.HiFi.sh
## hifiasm.sh
``` bash
hifiasm -o PigNoseTurtle -t 60 zb_1.fq  zb_2.fq
```
#02.Hi-C.sh
## mapping_arima.sh
``` bash
#! /bin/bash

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

# main path
main_folder="/public/home/zhengjiangmin/Software/mapping_pipeline"

SRA='2500_raw'
LABEL='zb'
BWA='/public/home/zhengjiangmin/Software/bwa/bwa'
SAMTOOLS='/public/home/zhengjiangmin/Software/samtools/bin/samtools'
IN_DIR='fastq'
REF='zb.fa'
FAIDX='$REF.fai'
PREFIX='zb'
RAW_DIR='out/bams'
FILT_DIR='out/filtered/bams'
FILTER=${main_folder}'/filter_five_end.pl'
COMBINER=${main_folder}'/two_read_bam_combiner.pl'
STATS=${main_folder}'/get_stats.pl'
PICARD='/public/home/zhengjiangmin/Software/picard/picard.jar'
TMP_DIR='/mnt/sata12'
PAIR_DIR='out/paired/bams'
# REP_DIR='/path/to/where/you/want/deduplicated/files'
REP_LABEL=$LABEL\_rep1
# MERGE_DIR='/path/to/final/merged/alignments/from/any/biological/replicates'
MAPQ_FILTER=10
CPU=12

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
# [ -d $REP_DIR ] || mkdir -p $REP_DIR
# [ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
[[ -f ${PREFIX}.bwt ]] || $BWA index -a bwtsw -p $PREFIX $REF

echo "### Step 0.5: Building fai file"
[[ -f ${PREFIX}.fa.fai ]] || $SAMTOOLS faidx ${PREFIX}.fa

echo "### Step 1.A: FASTQ to BAM (1st)"
$BWA mem -t $CPU `basename -s .fa $REF` $IN_DIR/$SRA\_1.fastq | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
$BWA mem -t $CPU `basename -s .fa $REF` $IN_DIR/$SRA\_2.fastq | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/$SRA\_1.bam $FILT_DIR/$SRA\_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

# ###############################################################################################################################################################
# ###                                           How to Accommodate Technical Replicates                                                                       ###
# ### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                                                    ###
# ### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
# ### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
# ### The code below is an example of how to merge technical replicates.                                                                                      ###
# ###############################################################################################################################################################
# #	REP_NUM=X #number of the technical replicate set e.g. 1
# #	REP_LABEL=$LABEL\_rep$REP_NUM
# #	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as technical replicates
# #   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
# #	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

# echo "### Step 4: Mark duplicates"
# java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

# $SAMTOOLS index $REP_DIR/$REP_LABEL.bam

# perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

# echo "Finished Mapping Pipeline through Duplicate Removal"

# #########################################################################################################################################
# ###                                       How to Accommodate Biological Replicates                                                    ###
# ### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                              ###
# ### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
# ### The code below is an example of how to merge biological replicates.                                                               ###
# #########################################################################################################################################
# #
# #	INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as biological replicates
# #   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
# #
# #	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
# #
# #	$SAMTOOLS index $MERGE_DIR/$LABEL.bam

# # perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# # echo "Finished Mapping Pipeline through merging Biological Replicates"
```
## yahs.sh
``` bash
Species="zb"

[[ -f ${Species}.fa.fai ]] || samtools faidx ${Species}.fa

/public/home/zhengjiangmin/Software/yahs-1.1/yahs ${Species}.fa out/paired/bams/*.bam

(~/Software/yahs-1.1/juicer pre yahs.out.bin yahs.out_scaffolds_final.agp $Species.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=10 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

[[ -f yahs.out_scaffolds_final.sizes ]] || bioawk -c fastx '{print $name, length($seq)}' yahs.out_scaffolds_final.fa > yahs.out_scaffolds_final.sizes
(java -jar -Xmx32G ~/Software/juicer-1.6/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt out.hic.part yahs.out_scaffolds_final.sizes) && (mv out.hic.part out.hic)
```
#03.BUSCO.sh
##busco5_genome.sh
``` bash
genome=$1
kind=geno
database=sauropsida_odb10
core=15
sp=Carettochelys insculpta
mkdir -p $genome\.busco && cd $genome\.busco
busco -m $kind -i ../$genome -o busco -l $database -c $core --config config.ini
```
##busco5_genome.sh
``` bash
genome=$1
kind=geno
database=vertebrata_odb10
core=15
sp=Carettochelys insculpta
mkdir -p $genome\.busco && cd $genome\.busco
busco -m $kind -i ../$genome -o busco -l $database -c $core --config config.ini
```
##busco5_genome.sh
``` bash
genome=$1
kind=geno
database=tetrapoda_odb10
core=15
sp=Carettochelys insculpta
mkdir -p $genome\.busco && cd $genome\.busco
busco -m $kind -i ../$genome -o busco -l $database -c $core --config config.ini 
```
#04.last.sh
##lastdb.sh
``` bash
lastdb -P30 -c -u NEAR Rafetus_swinhoei.NEAR Rafetus_swinhoei.fa
```
##lastal.sh
``` bash
lastal -P20 -i2G -m10  pigNoseHomo.fa > pigNoseHomo.maf ;tail -n1 lastal_outdir/pigNoseHomo/pigNoseHomo.maf last-split -m10 pigNoseHomo.maf >pigNoseHomo.1.maf ;last-split -r -m10 pigNoseHomo.1.maf > pigNoseHomo.2.maf ;maf-sort pigNoseHomo.2.maf >pigNoseHomo.sort.maf
```
#05.Repeat_annotation.sh
##TRF.sh
``` bash
trf pigNoseTurtle.fa 2 7 7 80 10 50 500 -d -h -ngs > pigNoseTurtle.fa.dat; perl convertTRF2gff.pl pigNoseTurtle.fa.dat pigNoseTurtle.fa.gff
```
##Repeatmasker.sh
``` bash
RepeatMasker -pa 12 -species 'turtles' -nolow -norna -no_is -gff pigNoseTurtle.fa -dir pigNoseHomo
```
##repeatProteinMask.sh
``` bash
RepeatProteinMask -engine ncbi -noLowSimple -pvalue 1e-04 input.fa
```
##RepeatModeler.sh
``` bash
BuildDatabase -name pigNoseTurtle -engine ncbi pigNoseTurtle.fa 
RepeatModeler -database pigNoseTurtle -pa 10
RepeatMasker -e ncbi -lib consensi.fa.classified -gff pigNoseTurtle.fa -dir pigNoseTurtle
```
#06.Gene_annotation.sh
##Assembly_transcriptome.sh
``` bash
rnaspades.py -1 muscle_1.fq.gz -2 muscle_2.fq.gz -o muscle;
rnaspades.py -1 liver_1.fq.gz -2 liver_2.fq.gz -o liver;
rnaspades.py -1 kidney_1.fq.gz -2 kidney_2.fq.gz -o kidney;
rnaspades.py -1 spleen_1.fq.gz -2 spleen_2.fq.gz -o spleen;
rnaspades.py -1 lung_1.fq.gz -2 lung_2.fq.gz -o lung;
rnaspades.py -1 trachea_1.fq.gz -2 trachea_2.fq.gz -o trachea;
cat muscle_trans.fa liver_trans.fa kidney_trans.fa spleen_trans.fa lung_trans.fa trachea_trans.fa > transcripts.fasta
```
##elimination_redundancy.sh
``` bash
cd-hit-est -i transcripts.fasta -o transcripts.fasta -T 30 -M 1600 > cdhit.out
#predict_cds.sh
TransDecoder.LongOrfs -t cdhit.out;
TransDecoder.Predict -t cdhit.out
```
##denovo.sh
``` bash
augustus --species=human --uniqueGeneId=true  scaffolds/ctg000000.fa >ctg000000.gff
##Genwise
#blat.sh
blat pig_nosedTurtle.fa Chicken.pep -makeOoc=pig_nosedTurtle.fa.dnax.ooc -t=dnax -q=prot out.psl
blat pig_nosedTurtle.fa query.pep -ooc=pig_nosedTurtle.fa.dnax.ooc -maxIntron=5000000 -noHead -t=dnax -q=prot pig_nosedTurtle.psl;
```
##genwise.sh
``` bash
genewise query.fa -nosplice_gtag -trev pig_nosedTurtle.fa -silent -pretty -pseudo -gff -cdna -trans > gene.genewise
```
##EVM.pl
``` perl
perl partition_EVM_inputs.pl --genome ../data/genome.sm.fa --gene_predictions ../data/denovo.gff --protein_alignments ../data/homolog.gff --segmentSize 5000000 --overlapSize 10000 --partition_listing partitions_list.out 

perl write_EVM_commands.pl \
 --genome data/genome.sm.fa \
 --weighqts weights.txt \
 --gene_predictions data/denovo.gff \
 --protein_alignments data/homolog.gff \
 --output_file_name evm.out \
 --partitions partition/partitions_list.out > work.sh

perl convert_EVM_outputs_to_GFF3.pl \
	--partition partition/partitions_list.out \
	--output_file_name evm.out \
	--genome data/genome.sm.fa
#The weight configurations were as follows: the weight of Augustus software predicted by denovo was 2, and the weight of homologous notes and transcripts of different tissues was 10.
```
#07.Function_annotation.sh
##cog
``` bash
blastp -query pigNoseTurtle.fa -out pigNoseTurtle.pep.fa.cog.blast -num_alignments 10 -db cog_clean.fa -outfmt 0 -evalue 1e-5  -num_threads 20
```
##interpro
``` bash
interpro/interproscan-5.44-79.0/interproscan.sh -appl Gene3D,Pfam,SMART,SUPERFAMILY,PRINTS,Coils -td pigNoseTurtle.pep.temp -b pigNoseTurtle.pep.fa  -i pigNoseTurtle.pep.fa -iprlookup -goterms -f tsv -t p -dp
```
##kegg
``` bash
blastp  -query pigNoseTurtle.pep.fa -out pigNoseTurtle.pep.fa.kegg.blast -num_alignments 10 -db kegg_all_clean.fa  -outfmt 0 -evalue 1e-5-num_threads 20
```
##nr
``` bash
diamond blastp -q pigNoseTurtle.pep.fa -o pigNoseTurtle.pep.fa.nr.blast -k 10 -d NR.fasta -f 0 -e 1e-5 -p 20 -M 20G --quiet
```
##swissprot
``` bash
blastp  -query pigNoseTurtle.pep.fa -out pigNoseTurtle.pep.fa.swissprot.blast -num_alignments 10 -db uniprot_sprot.fasta -outfmt 0 -evalue 1e-5  -num_threads 20
```
##trembl
``` bash
diamond blastp -d uniprot_trembl.fasta -q pigNoseTurtle.pep.fa -o pigNoseTurtle.pep.fa.trembl.blast -k 10 -f 0 -e 1e-5 -p 20 -M 20G --quiet
```
````
