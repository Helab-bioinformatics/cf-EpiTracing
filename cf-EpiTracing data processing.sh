for i in `ls *_2.fq.gz`
do 
base=$(basename $i "_2.fq.gz") 
cutadapt -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B CTGTCTCTTATACACATCTGACGCTGCCGACGA -q 20 -O 10  --trim-n  -m 30  --max-n 0.1  -o  ./${base}_trimmed.R1.fq.gz  -p ./${base}_trimmed.R2.fq.gz  ${base}_1.fq.gz ${base}_2.fq.gz  
done

for sample in `ls *_trimmed.R1.fq.gz`; do base=$(basename $sample "_trimmed.R1.fq.gz");  bowtie2  -p 20 -X 2000 -x /media/helab/genomes/hg19_bowtie2/hg19 -1 ./${base}_trimmed.R1.fq.gz -2 ./${base}_trimmed.R2.fq.gz -S ./${base}.hg19.sam  2>./${base}.hg19.align.log; done
for sample in `ls *_trimmed.R1.fq.gz`; do base=$(basename $sample "_trimmed.R1.fq.gz");  bowtie2  -p 20 -X 2000 -x /media/helab/genomes/dm3_bowtie2/genome -1 ./${base}_trimmed.R1.fq.gz -2 ./${base}_trimmed.R2.fq.gz -S ./${base}.dm3.sam  2>./${base}.dm3.align.log; done

for sample in `ls *.sam`
do
base=$(basename $sample ".sam")
samtools view -h -b -S  ${base}.sam   |  samtools sort  -  |  samtools view -h -bq 1 -  >  ./${base}_sort_uni.bam
done

for i in `ls *_sort_uni.bam`
do 
base=$(basename $i "_sort_uni.bam")
java -jar /media/helab/software/picard-tools-2.2.4/picard.jar MarkDuplicates REMOVE_DUPLICATES=true   I=${base}_sort_uni.bam   O=./${base}_rmdup_picard.bam  M=./${base}_rmdup_picard.txt 
done

for i in *_rmdup_picard.bam
do
samtools index $i
done

for i in ./*_rmdup_picard.bam
do
base=$(basename $i "_rmdup_picard.bam")
num1=10000000
num2="$(samtools view -c  $i  2>&1 )"
res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
bamCoverage --scaleFactor  $res -b  $i  -e 300  --smoothLength 500 -p 20 -o  ./${base}.ext500.smo200.bw
done

for i in ./*_rmdup_picard.bam
do 
base=$(basename $i "_rmdup_picard.bam")
java -jar  /media/helab/software/picard-tools-2.2.4/picard.jar  CollectInsertSizeMetrics   I=$i    O=./${base}_insertMetrics.txt  H=./${base}_insertHistogram.pdf   M=0.5
done

for x in ./*_sort_uni.bam
do 
base=$(basename $x "_sort_uni.bam")
#samtools index $x
samtools view -c  ${base}_sort_uni.bam  >./${base}.uni.totolreads.log
done

for j in ./*_rmdup_picard.bam
do 
base=$(basename $j "_rmdup_picard.bam")
samtools view -c ${base}_rmdup_picard.bam  >./${base}.totolreads.log
done
