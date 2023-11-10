#remove_n_bases.py

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
