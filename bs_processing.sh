# Analysis of BS converted data

# Things to look at before starting
# 1. Paired end of single end?
# 2. RRBS or WGBS?
# 3. Special library prepartion? Single cell?
#  ^ https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html  (section 8)
# 4. What controls were used? Usually lambda.  
# 5. Do we have C -> T converted reference genome already?


############################################################
#                        Trimming                          #
############################################################
# Trim files of low quality reads and adapters
#$ -S /bin/bash
#$ -pe serial 1
#$ -l mfree=10G -l h_rt=100:00:00

# List samplenames no commas or quotes
ids=(
A
B
C
)

# Go to location of files
cd /net/hawkins/vol/OvCSC_2017/REST_OF_LOCATION

# Loop through files to trim
# trim galore options to be aware of:
# --fastqc : quality control, makes html report
# -q : quality cut off, worth lowering if mapping bad
#		Use 30 for BS
# --paired : takes 2 fastq files, return matching read pairs
# --gzip : output is compressed

for i in ${ids[@]}; do
	trim_galore --fastqc --paired -q 30  ${i}.fq
done

# Check trimming by looking at .sh.o file made from running script
# Check that an adapter was detected
# Expect more to be trimmed out that ChIP or RNA-Seq


############################################################
#                        Mapping                           #
############################################################
# Mapping BS data can take a long time
# First script will be run for each sample
# Second script will submit first script on each sample

#1 bismark_map.sh

#$ -S /bin/bash
#$ -pe serial 1
#$ -l mfree=10G -l h_rt=100:00:00

cd /net/hawkins/vol/OvCSC_2017/REST_OF_LOCATION

# Many options for bismark
# Worth talking to who made library
# --rrbs is important in removing 'fake' methylation!
# -L and -N specify length of substring and # of max mistmatches to align
#	20 and 0 are default
# -o is output directory
# long directory listed is location of bismark genome

#map to lambda
bismark --se --rrbs -L 20 -N 0 /net/hawkins/vol1/scripts_files/lambda_genome -o bismark_lambda_${file} ${file}_trimmed.fq

# map to hg19
bismark --se --rrbs -L 20 -N 0 /net/hawkins/vol1/scripts_files/Bisulfite_Genome_Bowtie2 -o bismark_hg19_${file} ${file}_trimmed.fq

# Check your .sh.o file for % methylation of lambda and hg19 
# lambda should be <5% in all contexts (CpG, CHG, CHH)
# Check hg19 for CpG ~75% and CHG and CHH being lower


#########################################################

#2    submit_scripts.sh

#$ -S /bin/bash
#$ -pe serial 1
#$ -l mfree=10G -l h_rt=100:00:00

cd /net/hawkins/vol/OvCSC_2017/REST_OF_LOCATION

ids=(
A
B
C
)

cd /net/hawkins/vol1/dux4_rnaseq/raw_fastq
for i in ${ids[@]}; do
  qsub -cwd -v file=${i} -N ${i} bismark_map.sh
done






