X-linked dystonia-parkinsonism (XDP) is a combined dystonia-parkinsonism syndrome, modified by a (CCCTCT)n repeat within the causal SINE-VNTR-Alu (SVA) retrotransposon insertion in the TAF1 gene. This GitHub repository provides bash scripts to estimate the hexanucleotide repeat length from long-read nanopore sequencing data of patients with XDP. Therefore, the software tools NCRF or STRique can be used.

1. Installation

git clone https://github.com/nanopol/xdp_bookchapter.git

2. Usage 

NCRF:

sbatch xdp_bookchapter/NCRF/ncrf.sh </path/to/sequencing/data/in_zip_format.zip> <PCR/CRISPR> <Identifier>

To execute the NCRF script, three input parameters are required i.e.,: the path to the zip folder containing the FASTQ and FAST5 files, the method you have used to enrich the target (PCR or CRISPR) and the identifier for the resulting output folder.

  
For example: sbatch xdp_bookchapter/NCRF/ncrf.sh /data/patient_1234.zip CRISPR XDP_1234

STRique: 

sbatch xdp_bookchapter/STRique/STRIQUE.sh </path/to/sequencing/data/in_zip_format.zip> <PCR/CRISPR> <Identifier>

To execute the STRique script, the three input parameters has to be provided analogous to NCRF.

  
For example: sbatch xdp_bookchapter/STRique/STRIQUE.sh /data/patient_1234.zip CRISPR XDP_1234

