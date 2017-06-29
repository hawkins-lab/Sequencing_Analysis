#####################################################################
#    Processing of Data for ChIP-Seq, MARS, anything with 'peaks'   #
#####################################################################


#####################################################################
#  						Demultiplexing					         	#
#####################################################################


#In terminal
cd /net/hawkins/vol1/NEXTSEQ_YYMMDD
#

# data are located in 
/net/fields/vol2/fieldslab-inst/nextseq/Output/<RUNFOLDER>/

# e.g
<RUNFOLDER> = 141016_NS500272_0021_AH1441BGXX <date>_<machine>_<run>_<FlowCellID>

# create SampleSheet.csv in run folder 

[Data],,
Sample_ID,Sample_Name,index
CV160,CV160,TGACCTTGA

# all samples must have same lenght barcode sequence

# create output folder "NEXTSEQ_YYMMDD" in /net/hawkins/vol1/

#$ -S /bin/bash
#$ -pe serial 1
#$ -l h_rt=100:00:00 -l mfree=5G


cd /net/hawkins/vol1/NEXTSEQ_141222
bcl2fastq -R /net/fields/vol2/fieldslab-inst/nextseq/Output/141222_NS500272_0031_AH3HKKBGXX/ -o /net/hawkins/vol1/NEXTSEQ_141222/ --barcode-mismatches 0 --no-lane-splitting

#####################################################################
#  					      	Processing					        	#
#####################################################################


# We will make 2 scripts

# 1. script with all commands for trimming, mapping, etc. 
#    use ${file} in place of all sample names

# 2. script will submit 1st script with ${file}

# So we don't overwrite our file
cp A.fastq.gz A_NAME.fastq.gz
cp B.fastq.gz B_NAME.fastq.gz

# basics of vi
# i insert mode
# ESC leave insert mode
# :wq save and exit
# :q! quit without saving
# dd delete line
# create file and open
vi script.sh

# make sure you are using a 

#####################################################################

# SCRIPT 1
#$ -S /bin/bash
#$ -pe serial 1
#$ -l h_rt=120:00:00 -l mfree=5G

 cd /net/hawkins/vol1/ChIP_walkthrough/raw_fq

# trim galore options to be aware of:
# --fastqc : quality control, makes html report
# -q : quality cut off, worth lowering if mapping bad
#		Use 30 for BS
# --paired : takes 2 fastq files, return matching read pairs
# --gzip : output is compressed
# paired end output will be _val_1.fq and _val_2.fq
trim_galore -q 30 --phred33 --gzip ${file}.fastq.gz 
gunzip ${file}_trimmed.fq.gz

# need to specify -1 and -2 for paired end
# -x : location of bowtie2 index
# -t: shows time taken for .o file
# -N: number mismatches in seed alignment
# -L: length of seed
bowtie2 -N 1 -L 25 -t -x /net/hawkins/vol1/scripts_files/hg19/hg19 -S ${file}.sam ${file}_trimmed.fq

# convert sam to binary alignment bam
samtools view -bS ${file}.sam > ${file}_nonSorted.bam

# sort bam file
samtools sort ${file}_nonSorted.bam ${file}

# create index for bam file
# some tools need this to run
samtools index ${file}.bam

# make bigwig file from bam
# bamCoverage is part of deepTools
# http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html
# consider ignoring duplicates,
# could there be duplicates not from amplification?
# ignore chrX since there could be 1 or 2!
bamCoverage --ignoreDuplicates --normalizeUsingRPKM --bam ${file}.bam --ignoreForNormalization chrX -o ${file}.bw

# save some space by compressing fastqs
# about 1/4 the size
gzip ${file}_trimmed.fq

# remove intermediate files
rm ${file}.sam
rm ${file}_nonSorted.bam
mv ${file}.bw /net/hawkins/vol1/ChIP_walkthrough/aligned
mv ${file}.bam* /net/hawkins/vol1/ChIP_walkthrough/aligned
mv ${file}_trimmed.fq.gz /net/hawkins/vol1/ChIP_walkthrough/trimmed_fq

####################################################################

# SCRIPT 2

#Submit script for each sample
#$ -S /bin/bash
#$ -pe serial 1
#$ -l h_rt=120:00:00 -l mfree=5G

# Change NAME to your name
indexnames=(
A_NAME
B_NAME
)

cd /net/hawkins/vol1/ChIP_walkthrough/
for i in ${indexnames[@]}; do
  qsub -cwd -v file=${i} -N ${i} chip_pipeline.sh
done

#####################################################################
#  					      	Visualization					      	#
#####################################################################

# Instruction for uploading are in new_ucsc_system.txt
