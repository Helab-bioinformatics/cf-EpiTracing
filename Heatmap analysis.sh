#Heatmap----TSS 5kb sort by gold standard
computeMatrix reference-point --referencePoint center -p 20 -R refGene_hg19_TSS5kb.bed  -S Gold_standard_raw.hg19.ext500.smo200.bw -a 5000 -b 5000 -bs 50 --sortRegions descend --missingDataAsZero -q -o Matrix.gz --outFileNameMatrix Matrix.tab --outFileSortedRegions Sort_region.bed 
nohup computeMatrix reference-point --referencePoint center -p 20 -R Sort_region.bed -S *.bw  -a 5000 -b 5000 -bs 50 --sortRegions keep --missingDataAsZero -q -o Matrix.bed.gz   --outFileNameMatrix  Matrix_all.tab   --outFileSortedRegions  Sort_region_2.bed 
nohup plotHeatmap -m Matrix.bed.gz -o Tss5kb.svg --sortRegions no --zMin 0  --zMax 5 --colorMap Purples --heatmapHeight 20  -x  all_merged  -y  TSS_up5kb-down5kb   --plotTitle    all
--removeOutliers
