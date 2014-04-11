##########################################################
# alignment with bwa-mem
##########################################################

export SET="SMS173 SMS193 SMS230 SMS231 SMS238 SMS239 SMS242 SMS243 SMS244 SMS253 SMS254 SMS255"
export DATA=/m/cphg-quinlan/cphg-quinlan/projects/sms-elsea/fastq
export OUT=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/bam
export STEPname=bmemSET
export BIN=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/bin
export GENOME=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/hg19/hg19_gatk.fa

echo "$BIN/bwa index -a bwtsw $GENOME" | qsub -q arq5xlab -W group_list=cphg_arq5x -V -l walltime=48:00:00 -l mem=20gb -l ncpus=8 -m e -M jmh2tt@virginia.edu > id

for sample in `echo $SET`
do
export QSUB="qsub -q arq5xlab -N $STEPname -W depend=afterany:$(cat id) -W group_list=cphg_arq5x -V -m be -M jmh2tt@virginia.edu -l mem=50gb -l ncpus=8";
echo "cd $OUT; $BIN/bwa mem -t 16 -M -R '@RG\tID:$sample\tSM:$sample' -v 1 $GENOME $DATA/$sample"_"1.fastq.gz $DATA/$sample"_"2.fastq.gz | $BIN/samtools view -q 20 -f 2 -Su - > $sample.bam" | $QSUB > $sample.id
done

###############################################################
# Sort bam
###############################################################

export SET="SMS173 SMS193 SMS230 SMS231 SMS238 SMS239 SMS242 SMS243 SMS244 SMS253 SMS254 SMS255"
export DIR=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/bam
export BIN=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/bin
export STEPNAME=sort
for sample in `echo $SET`
do
export QSUB="qsub -q arq5xlab -N $STEPNAME -W depend=afterany:$(cat $sample.id) -W group_list=cphg_arq5x -V -m be -M jmh2tt@virginia.edu -l mem=50gb -l ncpus=8";
echo "cd $DIR; $BIN/samtools sort -m 1G -@ 4 $sample.bam $sample.sorted" | $QSUB > $sample'2'.id
done


##################################################################
#bam index
##################################################################

export SET="SMS173 SMS193 SMS230 SMS231 SMS238 SMS239 SMS242 SMS243 SMS244 SMS253 SMS254 SMS255"
export OUTDIR=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/bam
export STEPNAME=bamindex
export BIN=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/bin
for sample in `echo $SET`;
do
export QSUB="qsub -q arq5xlab -W depend=afterany:$(cat $sample'2'.id) -W group_list=cphg_arq5x -V -N $STEPNAME -m bea -M jmh2tt@virginia.edu -l mem=50gb -l ncpus=8";
echo "cd $OUTDIR; $BIN/samtools index $sample.sorted.bam" | $QSUB > $sample'3'.id
done


###################################################################
# mark duplicates from sorted bam with PICARD
###################################################################

export SET="SMS173 SMS193 SMS230 SMS231 SMS238 SMS239 SMS242 SMS243 SMS244 SMS253 SMS254 SMS255"
export OUTDIR=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/markdups
export WORK=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/bam
export STEPNAME=dedup
export BIN=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/bin/picard-tools-1.105
for sample in `echo $SET`;
do
export QSUB="qsub -q arq5xlab -W depend=afterany:$(cat $sample'3'.id) -W group_list=cphg_arq5x -V -N $STEPNAME -m bea -M jmh2tt@virginia.edu -l mem=50gb -l ncpus=8";
echo "cd $OUTDIR; java -Xmx4g -jar $BIN/MarkDuplicates.jar \
INPUT=$WORK/$sample.sorted.bam \
OUTPUT=$sample.sorted.markdup.bam \
TMP_DIR=$OUTDIR \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true \
CREATE_INDEX=true \
METRICS_FILE=$sample.markdup_metrics" | $QSUB > $sample'4'.id
done

##########################################################
## Variant calling with Freebayes
##########################################################

export SMSHOME=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/markdups
export VCF=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe
export GENOME=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/hg19/hg19_gatk.fa
export SET="SMS173 SMS193 SMS230 SMS231 SMS238 SMS239 SMS242 SMS243 SMS244 SMS253 SMS254 SMS255"
export NAME=sms
cd $SMSHOME
echo $SET | tr ' ' '\n' | awk -v var=$SMSHOME '{print var"/"$0".sorted.markdup.bam"}' > bam_list
export STEPNAME=freebaysmsdedup
export BIN=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/bin
for chrom in {1..22} X Y
do
export QSUB="qsub -W group_list=cphg_arq5x -q arq5xlab -V -l select=1:mem=8gb:ncpus=1 -N $STEPNAME -m bea -M jmh2tt@virginia.edu";
echo "$BIN/freebayes -L $SMSHOME/bam_list -f $GENOME -0 --min-coverage 10 \
--region $chrom \
[[[****OTHER ARGUMENTS HERE***]]]]
> $VCF/$NAME.rmdup.freebayes.$chrom.vcf" | $QSUB
done

(cat freebayes.1.vcf; cat freebayes*.vcf | grep -v ^# | awk '$1 != "1"') > $VCF/$NAME.rmdup.freebayes.genome.vcf

#######################################################
## Annotate VCF with snpEff
#######################################################

export SNPEFF=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/software/snpEff.jar
export VCF=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/vcfs
export STEPNAME=annotate

echo "java -jar $SNPEFF -i vcf -o vcf GRCh37.75 $VCF/sms.rmdup.all.freebayes.vcf > $VCF/sms.annotated.vcf" | qsub -V -q arq5xlab -N $STEPNAME -W group_list=cphg_arq5x -l mem=50gb -l ncpus=8

#if the field length is a problem: cat sms.rmdup.all.freebayes.vcf | awk '{ if (length($4) == length($5)) print }' > sms.modified.vcf

#######################################################
## Annotate VCF with VEP
#######################################################

export VEP=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vep
export VCF=/m/cphg-quinlan3/cphg-quinlan3/jmh2tt/vcfpipe/vcfs
#export STEPNAME=annotate

perl $VEP/ensembl-tools-release-75/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $VCF/sms.rmdup.all.freebayes.vcf \
--dir_cache $VEP/.vep/ --sift b --polyphen b --symbol --numbers --total_length -o $VCF/output.vcf --vcf --cache