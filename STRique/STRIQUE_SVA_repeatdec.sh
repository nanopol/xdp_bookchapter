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
#SBATCH --mail-user=
#  Find your job easier with a name:
#SBATCH --job-name=NCRF_SVA

#Initialize the module system:
source /etc/profile.d/modules.sh
# Allow aliases (required by some modules):
shopt -s expand_aliases
# Load your necessary modules:
module load strique/0.4.2
module load minimap2

#Define the home directory and path to sequencing data, your config file and input folders: 
HOME_DIR=/home/yourname
CONFIG=$HOME_DIR/config_SVA_$2.tsv
INPUT=$SCRATCH/inputdata
SVA=$SCRATCH/SVA
RESULTS=$SCRATCH/results
FAST5=$SCRATCH/fast5_files
SEQ_DATA=/data/youdatafolder

mkdir $INPUT
mkdir $SVA
mkdir $RESULTS
mkdir $FAST5


#Copy sequencing data and reference file to input:
unzip $SEQ_DATA/$1.zip -d  $INPUT
cp -a $INPUT/*.fastq $SVA
cp -a $HOME_DIR/SVA_$2.fa $SVA

#Make file of file names with STRiqe:
cp -a $HOME_DIR/SVA_$2.fa $INPUT
python $STRIQUE/scripts/STRique.py index $INPUT > $INPUT/SVA_reads.fofn
cp -a $CONFIG $INPUT

#Make alignment: 
cat $SVA/*.fastq > $SVA/SVA.fastq
minimap2 -a -x map-ont $SVA/SVA_$2.fa $SVA/SVA.fastq > $SVA/SVA.sam
samtools view -S -b $SVA/SVA.sam > $SVA/SVA.bam
samtools sort $SVA/SVA.bam -o $SVA/SVA.sorted.bam
samtools index $SVA/SVA.sorted.bam

#Make alignment with reads >1kb alignment length:
mkdir $SCRATCH/tempo/
cp -a $SVA/SVA.sorted.bam $SCRATCH/tempo/
samtools view -h $SCRATCH/tempo/SVA.sorted.bam | perl -lane '$l = 0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l < 1000 or /^@/' | samtools view -bS - > $INPUT/SVA_lenfilt.sorted.bam
samtools index $INPUT/SVA_lenfilt.sorted.bam

#Run STRique on aligned reads:
samtools view $INPUT/SVA_lenfilt.sorted.bam | python $STRIQUE/scripts/STRique.py count $INPUT/SVA_reads.fofn \
$STRIQUE/models/r9_4_450bps.model $INPUT/config_SVA_$2.tsv > $HOME_DIR/XDP_$2_$1.tsv

