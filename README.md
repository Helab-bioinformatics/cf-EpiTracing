# cf-EpiTracing data processing.sh
This file provides data pre-treatment pipeline. Briefly, we filtered out sequencing reads by indexing barcode sequences with any mismatch. After removing adaptors and low-quality bases by cutadapt (v.1.11), paired-end cf-EpiTracing reads were mapped to the human reference genome hg19 and Drosophila reference genome dm3 using Bowtie2 (v.2.2.9). The uniquely mapped reads with map quality greater than 30 were used for the following analyses. PCR duplicates were removed by Picard (v.2.2.4) (http://broadinstitute.github.io/picard). Only uniquely mapped, non-duplicated reads were used for analysis. To visualize cf-EpiTracing signals, we used Deeptools (v.2.2.3) to normalize reads to one million and calculate coverage in continuous 50-bp bins in bam files and generated track files (bigwig format). 

# ChromHMM.sh
This file provides ChromHMM pipeline, including the building up of 18-state ChromHMM model and chromatin state annotation of plasma samples using histone modifications of cf-chromatin.

# Correlation analysis.sh
This file provides quality control pipeline for calculating correlation of tested samples with reference.

# Heatmap analysis.sh
This file provides quality control pipeline for visualizing heatmap similarity of tested samples with reference.

# Identifying tissue specific regions.R
This file provides pipeline of identifying tissue specific signatures for tissues and primary cells using ChromHMM annotataion in Rstudio.

# Combining chromatin states.R
This file provides pipeline of combining chromatin state annotations of tested samples in Rstudio.

## Bed files provides with genome reference data that are necessary in analyses.
