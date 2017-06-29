#####################################################################
#    Analysis for for ChIP-Seq, MARS, anything with 'peaks'         #
#####################################################################

#1. Peaks, regions, overlaps with bed files

#2. Plots/stats with alignments

#####################################################################
#1. Calling Peaks

qlogin -l mfree=10G
cd /net/hawkins/vol1/ChIP_walkthrough/aligned/

# Using mac14 peak caller find regions of enrichment
# -t treatment, what you ChIP'd
# -c control, input: 
# --nomodel for histone modifications
# --nolambda faster doesn't make local background
# -S outputs wig, don't need since we have bigwig
# -g genome size, has hs, mm built in
#	you would need a number for other genomes
macs14 -t EH123.bam -n test_macs -c A.bam --nomodel --nolambda -S -g hs
# bed file output from macs14 is your peaks
# can be upload to UCSC same as bigwig files

# The extra columns from macs14 on intensity can be removed
awk  '{print $1 "\t" $2 "\t" $3}' test_macs_peaks.bed > output_peaks.bed

# bedtools is a suite of tools for common operations 
# http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html

# Some common ones

#sort 
bedtools sort -i combo.bed > sorted.bed

# merge
# will "flatten" bed files
# needs to be presorted
bedtools merge -i sorted.bed > merged.bed

# slop
# get additional flanking regions
# -b bp on both sides
# -l bp from start
# -r bp from end
# -s will use + and - to work directionally 
# -g genome chromosome size file
# prevents going outside a chromosome's coordinates
bedtools slop -i merged.bed -b 100 -g /net/hawkins/vol1/scripts_files/hg19_chromInfo.txt > with_flanking.bed

# intersect
# -a file1
# -b file2, can be multiple files
# -wa return complete regions of A with overlap
# -u do not repeat regions overlapping >1 time
# -f minimum fraction for overlap as decimal

# just subsection of peak overlapping LMR
bedtools intersect -a  output_peaks.bed -b lmrs.bed > peak_lmr_subset.bed
# keep whole peak with any LMR overlap
bedtools intersect -a  output_peaks.bed -b lmrs.bed -u > peak_lmr_whole.bed
# only keep whole peak if >50% overlaps LMR
bedtools intersect -a  output_peaks.bed -b lmrs.bed -u -f 0.5 > peak_lmr_whole_0.5.bed

# Look at length of bed files
wc -l *.bed

#####################################################################
# 2. deepTools
# Useful tools for analysis and making figures
# http://deeptools.readthedocs.io/en/latest/content/list_of_tools.html

# go back to bam files
cd ..

# multiBamSummary
# Get read counts of bam file across bins or bed file regions
# useful for comparisons and quantification
# -bs bin size
multiBamSummary bins -bs 10000 --bamfiles A.bam B.bam -out readCounts.npz 
# change bins to BED-file to quantify across bed regions
multiBamSummary BED-file --BED test_macs/peak_lmr_whole.bed \
--bamfiles A.bam B.bam  EH123.bam -out readCounts2.npz \
--outRawCounts readcounts2.tsv 

# look at output for readcounts2.tsv 
head readcounts2.tsv 
# use these readcounts for analysis
# make sure to normalize to sample depth and width of region

# The .npz file can be used to make plots with deepTools
plotCorrelation -in readCounts2.npz  --corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation" \
--whatToPlot heatmap -o heatmap_AB.png

# Heatmaps from bigwigs
# Need bed file and at least 1 bigwig
# scale-regions makes all widths appear equal
# reference-point goes set distance from point
computeMatrix scale-regions --regionsFileName test_macs/peak_lmr_whole.bed  --scoreFileName A.bw B.bw /net/hawkins/vol1/CD8_Benaroya/EH123.bw --binSize 20 --missingDataAsZero --sortRegions descend --outFileName test_matrix.npz

# plotHeatmap has a ton of options
# http://deeptools.readthedocs.io/en/latest/content/tools/plotHeatmap.html
-m test_matrix.npz  -out test_HM.pdf  --colorMap YlGnBu  \
--whatToShow 'heatmap and colorbar' --kmeans 3










