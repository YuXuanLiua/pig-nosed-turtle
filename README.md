# pig-nosed-turtle
A chromosome-level genome assembly of the pig-nosed turtle (Carettochelys insculpta)
## soapec.sh
``` bash
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 21 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 23 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 25 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
/public/home/zhuchenglong/software/soapec/SOAPec_bin_v2.03/bin/KmerFreq_HA -k 27 -t 50 -p PigNoseTurtle -l reads.list -L 150 >soapec.kmerfreq.log 2>soapec.kmerfreq.err
```
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
##busco5_genome.sh_sauropsida_odb10
``` bash
genome=$1
kind=geno
database=/data01/liuyuxuan2023/project/Busco_database/sauropsida_odb10
core=15
sp=zebrafish
rm -rf $genome\.busco
mkdir -p $genome\.busco && cd $genome\.busco

source /public/software/profile.d/compiler_gnu-7.2.0.sh
#augustus
export AUGUSTUS_CONFIG_PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/config
export PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/bin:$PATH
export PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/scripts:$PATH

export PYTHONPATH=/public/home/bioworker/software/BUSCO/sepp-4.5.1.sepp/lib/python3.9/site-packages:/home/software/Python-3.9.6/build/lib/python3.9/site-packages:/public/home/zhuchenglong/software/busco/build-5.2.2/lib/python3.9/site-packages/:$PYTHONPATH
export PATH=/home/software/Python-3.9.6/build/bin:/public/home/zhuchenglong/software/busco/pplacer-Linux-v1.1.alpha19:/public/home/zhuchenglong/software/metaeuk/metaeuk-4-a0f584d/build/bin/:$PATH

echo Start at: `date`
/data01/liuyuxuan2023/project/Software/busco/busco-5.2.2/bin/busco -m $kind -i ../$genome -o busco -l $database -c $core --config /public/home/lichunhui/software/busco/busco-5.2.2/config/config.ini
echo End   at: `date`
```
##busco5_genome.sh_vertebrata_odb10
``` bash
genome=$1
kind=geno
database=/data01/liuyuxuan2023/project/Busco_database/vertebrata_odb10
core=15
sp=zebrafish
rm -rf $genome\.busco
mkdir -p $genome\.busco && cd $genome\.busco

source /public/software/profile.d/compiler_gnu-7.2.0.sh
#augustus
export AUGUSTUS_CONFIG_PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/config
export PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/bin:$PATH
export PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/scripts:$PATH

export PYTHONPATH=/public/home/bioworker/software/BUSCO/sepp-4.5.1.sepp/lib/python3.9/site-packages:/home/software/Python-3.9.6/build/lib/python3.9/site-packages:/public/home/zhuchenglong/software/busco/build-5.2.2/lib/python3.9/site-packages/:$PYTHONPATH
export PATH=/home/software/Python-3.9.6/build/bin:/public/home/zhuchenglong/software/busco/pplacer-Linux-v1.1.alpha19:/public/home/zhuchenglong/software/metaeuk/metaeuk-4-a0f584d/build/bin/:$PATH

echo Start at: `date`
/data01/liuyuxuan2023/project/Software/busco/busco-5.2.2/bin/busco -m $kind -i ../$genome -o busco -l $database -c $core --config /public/home/lichunhui/software/busco/busco-5.2.2/config/config.ini
echo End   at: `date`
```
##busco5_genome.sh_tetrapoda_odb10
``` bash
genome=$1
kind=geno
database=/data01/liuyuxuan2023/project/Busco_database/tetrapoda_odb10
core=15
sp=zebrafish
rm -rf $genome\.busco
mkdir -p $genome\.busco && cd $genome\.busco

source /public/software/profile.d/compiler_gnu-7.2.0.sh
#augustus
export AUGUSTUS_CONFIG_PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/config
export PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/bin:$PATH
export PATH=/public/home/lichunhui/software/augustus/augustus-3.3.3/scripts:$PATH

export PYTHONPATH=/public/home/bioworker/software/BUSCO/sepp-4.5.1.sepp/lib/python3.9/site-packages:/home/software/Python-3.9.6/build/lib/python3.9/site-packages:/public/home/zhuchenglong/software/busco/build-5.2.2/lib/python3.9/site-packages/:$PYTHONPATH
export PATH=/home/software/Python-3.9.6/build/bin:/public/home/zhuchenglong/software/busco/pplacer-Linux-v1.1.alpha19:/public/home/zhuchenglong/software/metaeuk/metaeuk-4-a0f584d/build/bin/:$PATH

echo Start at: `date`
/data01/liuyuxuan2023/project/Software/busco/busco-5.2.2/bin/busco -m $kind -i ../$genome -o busco -l $database -c $core --config /public/home/lichunhui/software/busco/busco-5.2.2/config/config.ini
echo End   at: `date`
```
##lastdb.sh
``` bash
/public/home/zhuchenglong/software/last/last-1469/bin/lastdb -P30 -c -u NEAR Rafetus_swinhoei.NEAR Rafetus_swinhoei.fa
```
##lastal.sh
``` bash
( /usr/bin/time -v /public/home/zhuchenglong/software/last/last-1282/bin/lastal -P20 -i2G -m10 -pHOXD70 Rafetus_swinhoei_genome_NEAR/Rafetus_swinhoei.NEAR genome_data/pigNoseHomo.fa > lastal_outdir/pigNoseHomo/pigNoseHomo.maf ) > lastal_outdir/pigNoseHomo/pigNoseHomo.maf.log 2>&1 ; tail -n1 lastal_outdir/pigNoseHomo/pigNoseHomo.maf >> lastal_outdir/pigNoseHomo/pigNoseHomo.maf.log ;/public/home/zhuchenglong/software/last/last-1282/bin/last-split -m10 lastal_outdir/pigNoseHomo/pigNoseHomo.maf > lastal_outdir/pigNoseHomo/pigNoseHomo.1.maf ;/public/home/zhuchenglong/software/last/last-1282/bin/last-split -r -m10 lastal_outdir/pigNoseHomo/pigNoseHomo.1.maf > lastal_outdir/pigNoseHomo/pigNoseHomo.2.maf ;/public/home/zhuchenglong/software/last/last-1282/bin/maf-sort lastal_outdir/pigNoseHomo/pigNoseHomo.2.maf > lastal_outdir/pigNoseHomo/pigNoseHomo.sort.maf
```
##TRF.pl.sh
``` bash
/public/home/fengchenguang/software/RepeatMasker/TRF/trf /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/../00.data/split_1000/ZhuBiGui.part-01.fa 2 7 7 80 10 50 500 -d -h -ngs > /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/trf_output/ZhuBiGui.part-01.fa.dat; /usr/bin/perl /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/scripts/convertTRF2gff.pl /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/trf_output/ZhuBiGui.part-01.fa.dat /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/trf_output/ZhuBiGui.part-01.fa.gff
```
##Repeatmasker.pl.sh
``` bash
/public/home/fengchenguang/software/RepeatMasker/TRF/trf /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/../00.data/split_1000/ZhuBiGui.part-01.fa 2 7 7 80 10 50 500 -d -h -ngs > /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/trf_output/ZhuBiGui.part-01.fa.dat; /usr/bin/perl /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/scripts/convertTRF2gff.pl /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/trf_output/ZhuBiGui.part-01.fa.dat /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/01.TRF/trf_output/ZhuBiGui.part-01.fa.gff
```
##repeatProteinMask.pl.sh
``` bash
cd /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/03.repeatProMasker/rmdir ;/public/home/fengchenguang/software/RepeatMasker/RepeatMasker/RepeatProteinMask -engine ncbi -noLowSimple -pvalue 1e-04/data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/03.repeatProMasker/repeatProteinMasker_output/ZhuBiGui.part-01/input.fa
```
##RepeatModeler.sh
``` bash
mkdir -p database ; cd database ; ln -s /data01/liuyuxuan2023/project/data/ZhuBiGui.fa ; /public/home/fengchenguang/software/RepeatModeler/RepeatModeler-open-1.0.11/BuildDatabase -name genome -engine ncbi genome.fa ; cd -
/public/home/fengchenguang/software/RepeatModeler/RepeatModeler-open-1.0.11/RepeatModeler -database ./database/genome -pa 10
cd /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/04.repeatModeler/rmdir; /public/home/fengchenguang/software/RepeatMasker/RepeatMasker/RepeatMasker -e ncbi -lib /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/04.repeatModeler/RM_database/consensi.fa.classified -gff /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/04.repeatModeler/../00.data/split_1000/ZhuBiGui.part-01.fa -dir /data01/liuyuxuan2023/project/templete/00.repeat/04.zhubigui/04.repeatModeler/repeatmasker_output/ZhuBiGui.part-01
```
##cog
``` bash
time /public/home/humingliang/software/blast/ncbi-blast-2.7.1+/bin/blastp -query /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/pigNoseHomo.pep.cut/pigNoseHomo.pep.1.fa -out /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/cog/pigNoseHomo.pep.1.fa.cog.blast -num_alignments 10 -db /public/home/humingliang/database/gene_function/cog/cog_clean.fa -outfmt 0 -evalue 1e-5  -num_threads 20
```
##interpro
``` bash
export PATH=/public/home/humingliang/software/java/jdk-11.0.7/bin:$PATH;export CLASSPATH=.:/public/home/humingliang/software/java/jdk-11.0.7/lib/dt.jar:/public/home/humingliang/software/java/jdk-11.0.7/lib/tools.jar;cd /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/interpro && time /public/home/humingliang/software/gene_function/interpro/interproscan-5.44-79.0/interproscan.sh -appl Gene3D,Pfam,SMART,SUPERFAMILY,PRINTS,Coils -td /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/interpro/pigNoseHomo.pep.1.temp -b /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/interpro/pigNoseHomo.pep.1.fa  -i /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/pigNoseHomo.pep.cut/pigNoseHomo.pep.1.fa -iprlookup -goterms -f tsv -t p -dp
```
##kegg
``` bash
time /public/home/humingliang/software/blast/ncbi-blast-2.7.1+/bin/blastp  -query /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/pigNoseHomo.pep.cut/pigNoseHomo.pep.1.fa -out /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/kegg/pigNoseHomo.pep.1.fa.kegg.blast -num_alignments 10 -db /public/home/humingliang/database/gene_function/kegg/kegg_all_clean.fa  -outfmt 0 -evalue 1e-5-num_threads 20
```
##nr
``` bash
time /public/home/humingliang/software/diamond/diamond-2.0.4/diamond blastp -q /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/pigNoseHomo.pep.cut/pigNoseHomo.pep.1.fa -o /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/nr/pigNoseHomo.pep.1.fa.nr.blast -k 10 -d /public/home/humingliang/database/gene_function/NR/NR_diamond/NR.fasta -f 0 -e 1e-5 -p 20 -M 20G --quiet
```
##swissprot
``` bash
/public/home/humingliang/software/blast/ncbi-blast-2.7.1+/bin/blastp  -query /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/pigNoseHomo.pep.cut/pigNoseHomo.pep.1.fa -out /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/swissprot/pigNoseHomo.pep.1.fa.swissprot.blast -num_alignments 10 -db /data01/liuyuxuan2023/project/templete/06.gene_function/uniprot_sprot.fasta -outfmt 0 -evalue 1e-5  -num_threads 20
```
##trembl
``` bash
time /public/home/humingliang/software/diamond/diamond-2.0.4/diamond blastp -d /public/home/humingliang/database/gene_function/trembl/uniprot_trembl.fasta -q /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/pigNoseHomo.pep.cut/pigNoseHomo.pep.1.fa -o /data01/liuyuxuan2023/project/templete/06.gene_function/03.ZhuBiGui/trembl/pigNoseHomo.pep.1.fa.trembl.blast -k 10 -f 0 -e 1e-5 -p 20 -M 20G --quiet
```





