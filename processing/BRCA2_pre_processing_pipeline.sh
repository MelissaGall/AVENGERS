#!/bin/bash
#Variables definition
seq_type="NS" #or 
seq_dir="/mnt/f/BRCA2_NTD/data"
out_dir="/mnt/f/BRCA2_NTD/qc"
#Run fastq
mkdir $out_dir
mkdir $out_dir/fastq
mkdir $out_dir/fastq/fastqc_out
fastqc $seq_dir/*.fastq.gz -o $out_dir/fastq/fastqc_out
echo "Running fastqc."
#Run seqprep 
cd $seq_dir
for f in *.fastq.gz; do mv "$f" "X$f"; done
rename 's/preRP/PRERP/' *.fastq.gz
rename 's/d3r/PRERP/' *.fastq.gz
rename 's/D3r/PRERP/' *.fastq.gz
rename 's/d3R/PRERP/' *.fastq.gz
rename 's/d3R/PRERP/' *.fastq.gz
rename 's/D3R/PRERP/' *.fastq.gz
rename 's/CISR/CISRP/' *.fastq.gz
rename 's/DMSOR/DMSORP/' *.fastq.gz
rename 's/OLAR/OLARP/' *.fastq.gz
rename 's/CisR/CISRP/' *.fastq.gz
rename 's/CIsR/CISRP/' *.fastq.gz
rename 's/OlaR/OLARP/' *.fastq.gz
mkdir -p Seqprep
mkdir -p Seqprep/R1
mkdir -p Seqprep/R2
mkdir -p Seqprep/merged
python3 /mnt/f/BRCA2_NTD/scripts/run_seqprep.py
echo "Running Seqprep on all samples in parallel."
sh run_seqprep.sh
#Remove n bases
cd Seqprep/merged
mkdir no_Ns
echo "Remove Ns"
python3 /mnt/f/BRCA2_NTD/scripts/run_remove_n_bases.py
sh run_remove_n_bases.sh
cd no_Ns
#Run needle
rename 's/^(X[0-9]*[a-z]*)/\1-/' *.fastq
mkdir sam
python3 /mnt/f/BRCA2_NTD/scripts/run_needle_to_sam.py
sh run_needle_to_sam.sh