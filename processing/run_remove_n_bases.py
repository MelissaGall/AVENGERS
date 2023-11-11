#run_remove_n_bases.py
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
shell_file = open('run_remove_n_bases.sh', 'w')
shell_file.write('cd ' + working_dir + '\n')

for i in os.listdir(working_dir):
    if i.endswith(".fastq.gz"):
        index_of_first_dot = i.find('.')
        index_of_first_under = i.find('_')
        index_of_extension = i.find('.fastq.gz')
        sample_name = i[:index_of_first_dot]
        sample_lib = i[:index_of_first_under]
        before_extension_name = i[:index_of_extension]
        shell_file.write("python3 /mnt/f/BRCA2_NTD/scripts/remove_n_bases.py " + i + ' ' + 'no_Ns/' + before_extension_name + '.fastq &\n')
    else:
        pass

shell_file.write('wait\n')
shell_file.close()