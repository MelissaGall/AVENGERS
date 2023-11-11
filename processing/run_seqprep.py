#run_seqprep.py
#remove_n_bases.py
#Original script from: 
#Findlay GM, Daza RM, Martin B, Zhang MD, Leith AP, Gasperini M, Janizek JD, Huang X, Starita LM, Shendure J. 
#Accurate classification of BRCA1 variants with saturation genome editing. 
#Nature. 2018 Oct;562(7726):217-222. 
#doi: 10.1038/s41586-018-0461-z. Epub 2018 Sep 12. PMID: 30209399; PMCID: PMC6181777.
#https://github.com/shendurelab/saturationGenomeEditing_pipeline
import os
import sys
import subprocess

working_dir = os.getcwd()
print(working_dir)
command_file_name = working_dir+"/run_seqprep.sh"
command_file = open(command_file_name, 'w')

for i in os.listdir(os.getcwd()):
    if i.endswith("R1_001.fastq.gz"): 
        index_of_first_under = i.find('_')
        index_of_R1 = i.find('R1')
        sample_name = i[:index_of_first_under]
        file_name = i[:index_of_R1]
        print(file_name+r'R1_001.fastq.gz')
        command_file.write(r'seqprep -f '+ file_name + r'R1_001.fastq.gz -r '+ file_name + r'R2_001.fastq.gz -1 Seqprep/R1/' + sample_name + r'.R1.fastq.gz -2 Seqprep/R2/' + sample_name + r'.R2.fastq.gz -A GGTTTGGAGCGAGATTGATAAAGT -B CTGAGCTCTCTCACAGCCATTTAG -M 0.1 -s Seqprep/merged/' + sample_name + '.merged.fastq.gz -m 0.001 -q 20 -o 20 &\n')
    else:
        pass
command_file.write("wait\n")
command_file.close()

#notes for seqprep
'''
adapter trimming for A and B param's:
-A GGTTTGGAGCGAGATTGATAAAGT #if PU1R is used (as seen in R1):  
-B CTGAGCTCTCTCACAGCCATTTAG #if PU1L is used (as seen in R2): 
-M 0.1
'''
