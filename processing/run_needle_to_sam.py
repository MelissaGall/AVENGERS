#run_needle_to_sam_BRCA1_pipeline.py
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
shell_file = open('run_needle_to_sam.sh', 'w')

for i in os.listdir(working_dir):
    print(i)
    if i.endswith(".fastq"): 
        index_of_first_dot = i.find('.')
        sample_name = i[:index_of_first_dot]
        index_of_first_dash = i.find('-')
        if 'r' in sample_name: #for i.e. X17r1-pre
            sample_amplicon = sample_name[:i.find('r')]
            print(sample_amplicon)
        else:
            sample_amplicon = i[:index_of_first_dash]
            print(sample_amplicon)
        sample_ref = '/mnt/f/BRCA2_NTD/reference/'+sample_amplicon+'.fa'
        shell_file.write("needleall -asequence " + sample_ref+ " -bsequence "+i+ " -gapopen 10 -gapextend 0.5 -outfile sam/" + sample_name[:index_of_first_dash]+'_'+sample_name[index_of_first_dash+1:] +".sam -aformat sam &\n")       
    else:
        pass
shell_file.write("wait\n")
shell_file.close()