#! /bin/bash
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=shortterm
#  Use one node:
#SBATCH --nodes=1
#  Request 10 cores (hard constraint):
#SBATCH -c 20
#  Request 10GB of memory (hard constraint):
#SBATCH --mem=100GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=1-0:0:0
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=100G
#  Notify me at job start and end:
#SBATCH --mail-type=ALL
#  Send the notifications to:
#SBATCH --mail-user=theresa.lueth@gmx.de
#  Find your job easier with a name:
#SBATCH --job-name=NCRF_SVA

#Initialize the module system:
source /etc/profile.d/modules.sh
# Allow aliases (required by some modules):
shopt -s expand_aliases
# Load your necessary modules:
module load ncrf/v1.01.02
module load minimap2
module load samtools

#Define the repeat motif, NCRF input folder, home directory and path to sequencing data in .zip format: 
IDENTIFIER=$3
MOTIF=AGAGGG
CONFIG=$HOME/ncrf_results_$3_$2
INPUT=$SCRATCH/inputdata
SEQ_DATA=$1
mkdir $INPUT
mkdir $CONFIG

#Copy sequencing data and reference file to input:
cp -a $HOME/xdp_bookchapter/SVA_$2.fa $INPUT
unzip $1 -d  $INPUT

#Prepare a BAM file containing only reads with >1kb alignment length:
cat $INPUT/*.fastq > $INPUT/$1.fastq
minimap2 -a -x map-ont $INPUT/SVA_$2.fa $INPUT/$3.fastq > $INPUT/$3.sam
samtools view -S -b $INPUT/$3.sam > $INPUT/$3.bam
samtools sort $INPUT/$3.bam -o $INPUT/$3.sorted.bam
samtools index $INPUT/$3.sorted.bam
samtools view -F 2308  $INPUT/$3.sorted.bam -bS > $INPUT/aligned_$3.sorted.bam
samtools view -h $INPUT/aligned_$3.sorted.bam | perl -lane '$l = 0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l > 1000 or /^@/' | samtools view -bS - > $INPUT/aligned_length_filt_$3.sorted.bam
samtools index $INPUT/aligned_length_filt_$3.sorted.bam

#Copy generated BAM and index file to results folder:

cp -a $INPUT/aligned_length_filt_$3.sorted.bam $CONFIG
cp -a $INPUT/aligned_length_filt_$3.sorted.bam.bai $CONFIG

#Convert the BAM file first to FASTQ and then to FASTA:
samtools fastq $INPUT/aligned_length_filt_$3.sorted.bam > $INPUT/aligned_length_filt_$3.fastq
cat $INPUT/aligned_length_filt_$3.fastq | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $INPUT/aligned_length_filt_$3.fasta


#Do NCRF analysis on FASTA file with set max. noise and min. length parameters:
cat $INPUT/aligned_length_filt_$3.fasta \
| NCRF $MOTIF --stats=events --positionalevents --maxnoise=20%  --minlength=108 \
| ncrf_sort.py --sortby=mratio  \
| tee $CONFIG/"${IDENTIFIER}_${MOTIF}_raw.summary" \
| ncrf_consensus_filter.py \
| ncrf_sort.py --sortby=mratio  \
| tee $CONFIG/"${IDENTIFIER}_${MOTIF}_refined.summary" \
| ncrf_summary.py \
> $CONFIG/"${IDENTIFIER}_${MOTIF}_summary.summary"

echo "--->NCRF analysis is finished"








