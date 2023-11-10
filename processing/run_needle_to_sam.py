#run_needle_to_sam_BRCA1_pipeline.py

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