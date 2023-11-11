#remove_n_bases.py
#Original script from: 
#Findlay GM, Daza RM, Martin B, Zhang MD, Leith AP, Gasperini M, Janizek JD, Huang X, Starita LM, Shendure J. 
#Accurate classification of BRCA1 variants with saturation genome editing. 
#Nature. 2018 Oct;562(7726):217-222. 
#doi: 10.1038/s41586-018-0461-z. Epub 2018 Sep 12. PMID: 30209399; PMCID: PMC6181777.
#https://github.com/shendurelab/saturationGenomeEditing_pipeline
import gzip
import sys
from itertools import islice
out_file = open(sys.argv[2],'w')

with gzip.open(sys.argv[1], 'rt') as my_reads:
    reads = 0
    reads_w_n = 0
    while True:
        my_fastq = list(islice(my_reads, 4))
        if not my_fastq:
            break
        else:
            reads+=1
            seq = my_fastq[1]
            if "N" in seq:
                reads_w_n += 1
            elif 'N' not in seq:
                out_file.write(''.join(my_fastq))
    print(sys.argv[1], 'finished processing.')
    print(reads_w_n, 'reads with N base')
    print(reads, 'total reads')
