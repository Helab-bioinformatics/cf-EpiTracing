#TSS 5-kb correlation
multiBigwigSummary BED-file -b *bw  -out Tss5kb.results.npz --BED   refGene_hg19_TSS5kb.bed --outRawCounts ReadCounts_correlation_tss5bkb.tab  -p  20 
plotCorrelation -in Tss5kb.results.npz  -o  Tss5kbcorrelation_heatmap.svg -c pearson --removeOutliers -p  heatmap  --skipZeros  --colorMap inferno
plotCorrelation -in Tss5kb.results.npz  -o  Tss5kbcorrelation_scatter.svg -c pearson --removeOutliers -p scatterplot  --skipZeros  --colorMap inferno

#Genome 10-kb correlation
multiBigwigSummary BED-file -b *bw  -out Genome10kb.results.npz --BED   hg19.10K.windows.bed --outRawCounts readCounts_correlation_genome10kb.tab  -p  20 
plotCorrelation -in Genome10kb.results.npz  -o  Genome10kb_correlation_heatmap.svg -c pearson --removeOutliers -p  heatmap  --skipZeros  --colorMap inferno
plotCorrelation -in Genome10kb.results.npz  -o  Genome10kb_correlation_scatter.svg -c pearson --removeOutliers -p scatterplot  --skipZeros  --colorMap inferno
