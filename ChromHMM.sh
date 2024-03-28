#Learn 18-state ChromHMM model with ChIP-seq data for tissues and primary cells
java -mx4000M -jar ChromHMM.jar BinarizeBed ./ChromHMM/CHROMSIZES/hg19.txt ChIPseq_reference_data.txt ./BinarizeBed_Reference
nohup java -mx40000M -jar ./ChromHMM/ChromHMM.jar LearnModel -p 18 ./BinarizeBed_reference ./Model_18_chromatin_state hg19

#Run Binary files for plasma samples (in 200-bp resolution)
java -mx4000M -jar ChromHMM.jar BinarizeBed ./ChromHMM/CHROMSIZES/hg19.txt Samplelist.txt ./BinarizeBed&

#Make annotation for plasma samples (in 200-bp resolution)
nohup java -mx100000M -jar ./ChromHMM/ChromHMM.jar MakeSegmentation -printposterior Model_18_chromatin_state.txt ./ChromHMM/BinarizeBed ./MakeSegment&

#Run Binary files for plasma samples (in 5-kb resolution)
java -mx4000M -jar ChromHMM.jar BinarizeBed -b 5000 ./ChromHMM/CHROMSIZES/hg19.txt Samplelist.txt ./BinarizeBed_5kb&

#Make annotation for plasma samples (in 5-kb resolution)
nohup java -mx100000M -jar ./ChromHMM/ChromHMM.jar MakeSegmentation -printposterior Model_18_chromatin_state.txt ./ChromHMM/BinarizeBed_5kb ./MakeSegment_5kb&

#Segment annotation to 200-bp
for i in ./*.bed
do
base=$(basename $i ".bed")
bedtools makewindows -b $i -w 200 -i src > ${base}_200.txt
done

#Sort segments to the same order
for i in ./*_200.txt
do
base=$(basename $i "_200.txt")
sort -k1,1 -k2,2n $i > ${base}_sort.txt
done

