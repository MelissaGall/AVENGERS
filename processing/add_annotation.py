#!/usr/bin/env python
#coding: utf-8
#Script used to add ClinVar and c/g nomemclature to the count files
#Load libraries
import pandas as pd
import numpy as np
import os
import subprocess
import operator
import re
#Load files and get list of amplicons
all_files = os.listdir('./variant_counts_no_thresh')
txt_list = [x for x in all_files if '.txt' in x]
amplicon_list = [s.split('RP')[0] for s in txt_list]
amplicon_list = list(np.unique(amplicon_list))
print(amplicon_list)
# Create a dictionary with the reference sequences for each exon
fasta_dir = '../../../../../../BRCA2/reference'
ref_seqs = {}
for amplicon in amplicon_list:
    ref_file = open(fasta_dir+'/'+amplicon+'.fa', 'r')
    ref_header = ref_file.readline()
    ref_seq = ref_file.readline().strip()
    ref_seqs[amplicon] = ref_seq
    ref_file.close()
#Create a dictionary in which each amplicon in the file is a key that points to a list containing the info for that amplicon [5' HDR, 3' HDR, cut_pos, cut_end, mut_start, mut_end] 
edits_per_amp = {}
editing_info = open('../../../../../reference/BRCA2_editing_data.txt', 'r')
editing_info_header = editing_info.readline()
for line in editing_info:
    edits = line.strip().split()
    edits_per_amp[edits[0]] = edits[1:]
#Load ClinVar data for BRCA2
BRCA_data = pd.read_csv('../../../../..//reference/scores/clinvar_result_BRCA2.txt', sep = '\t', header = 0)
BRCA_data = BRCA_data[BRCA_data['Gene(s)'] == 'BRCA2']
#Load splicind data
splicing = pd.read_excel('../../../../../reference/splicing/all_exons_summary_spliceAI.xlsx')
#AA dictionary
DNA_Codons = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TGT': 'C', 'TGC': 'C',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'TTT': 'F', 'TTC': 'F',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAT': 'H', 'CAC': 'H',
    'ATA': 'I', 'ATT': 'I', 'ATC': 'I',
    'AAA': 'K', 'AAG': 'K',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATG': 'M',
    'AAT': 'N', 'AAC': 'N',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGA': '*'
}
#Define dictionary with the values used to convert g to c nomenclature (for GRCh38)
c_to_g_dict = {'2':  32316448,
               '3': 32319009,
              }
#Define function to get difference
def dif(a, b):
    return [i for i in range(len(a)) if a[i] != b[i]]
#Define function to find the difference in codon/AA
def seq_to_codons(seq):
    frame = (exon_start - genomic_pos) % 3
    seq = seq[frame:]
    prot_pos = prot_start - ((exon_start - genomic_pos) // 3) 
    ref_seq_var = list(ref_seq)
    ref_seq_var[int(edits3[0].split('-')[0]) - 1] = edits3[0].split('-')[2]
    ref_seq_var[int(edits5[0].split('-')[0]) - 1] = edits5[0].split('-')[2]
    ref_seq_var = ''.join(ref_seq_var)
    ref_seq_var = ref_seq_var[frame:]
    for i in range(0, len(seq), 3):
        codon_seq = seq[i:i+3]
        if len(codon_seq) < 3 or len(ref_seq_var[i:i+3]) < 3 :
            break
        seq_codon_ref_seq = ref_seq_var[i:i+3]
        if codon_seq != seq_codon_ref_seq:
            pos = dif(codon_seq, seq_codon_ref_seq)
            ref_nuc = seq_codon_ref_seq[pos[0]]
            var_nuc = codon_seq[pos[0]]
            ref_aa = DNA_Codons[seq_codon_ref_seq]
            var_aa = DNA_Codons[codon_seq]
            nuc_change = ref_nuc + ':' + var_nuc
            if prot_pos < prot_start or prot_pos > prot_end :
                prot_change = 'Intronic'
            else:
                prot_change = ref_aa + str(prot_pos) + var_aa
        prot_pos += 1
    return nuc_change, prot_change
#Run the script for each count file
for i in txt_list:
    #Load file and extract info on the amplicon
    df = pd.read_csv(('./variant_counts_no_thresh/' + i), sep = '\t', header = 0)
    file_str = i.split('.txt')[0]
    amplicon = i.split('RP')[0]
    exon = re.sub('[^0-9]', '', amplicon)
    c_to_g_val = c_to_g_dict[exon]
    ref_seq = ref_seqs[amplicon]
    my_editing_info = edits_per_amp[amplicon]
    genomic_pos = int(my_editing_info[6])
    prot_start = int(my_editing_info[7])
    exon_start = int(my_editing_info[8])
    exon_stop = int(my_editing_info[9])
    prot_end = int(my_editing_info[10])
    edits3 = my_editing_info[0].strip().split(',')
    edits5 = my_editing_info[1].strip().split(',')
    if my_editing_info[0] == my_editing_info[1] :
        PAM_seq = my_editing_info[0]
    else :
         PAM_seq = my_editing_info[0] + "," + my_editing_info[1]
    print("Start working on", file_str)
    #Remove WT and just PAM sequences
    df = df[df['edit_string'] != "-WT"]
    df = df[df['edit_string'] != PAM_seq]
    #Get position of the mutation and nucleotide change
    df['pos'] = df['edit_string'].str.replace(my_editing_info[0], '')
    df['pos'] = df['pos'].str.replace(my_editing_info[1], '')
    df['pos'] = df['pos'].str.replace(',', '')
    #df = df[df['pos'].str.len() >= 3]
    df['pos'] = df['pos'].str.partition('-')[0].astype(int)
    df['chr_pos'] = df['pos'] + genomic_pos - 1
    df['nuc_change'], df['AA_change'] = zip(*df['variant'].map(seq_to_codons))
    df['AA_change'] = np.where(df['chr_pos'] < exon_start, 'Intronic', df['AA_change'])
    df['AA_change'] = np.where(df['chr_pos'] > exon_stop, 'Intronic', df['AA_change'])
    #Load ClinVar data for the corresponding part of the gene to avoid working with huge data frame
    BRCA_data_subset = BRCA_data[BRCA_data.GRCh38Location >= genomic_pos]
    BRCA_data_subset = BRCA_data_subset[BRCA_data_subset.GRCh38Location < genomic_pos + 1000]
    BRCA_data_subset = BRCA_data_subset[['GRCh38Location', 'Protein change', 'dbSNP ID', 'Clinical significance (Last reviewed)', 'Canonical SPDI']]
    #Clean some data: nucleotide change, clinical significance...
    BRCA_data_subset['change'] = BRCA_data_subset['Canonical SPDI'].str[-3:]
    BRCA_data_subset['clinical_sign'] = BRCA_data_subset['Clinical significance (Last reviewed)'].str.split('(').str[0]
    BRCA_data_subset.GRCh38Location = BRCA_data_subset.GRCh38Location.astype('int')
    BRCA_data_subset = BRCA_data_subset[['GRCh38Location', 'dbSNP ID', 'change', 'clinical_sign', 'Protein change']]
    #Merge our data with ClinVar data
    final = df.merge(BRCA_data_subset, how = 'left', left_on = ['chr_pos', 'nuc_change'], right_on = ['GRCh38Location', 'change'])
    #Add c and g nomenclature
    final['c_val'] = final['chr_pos'] - c_to_g_val
    try:
        min_c_val = final.c_val[final.chr_pos == exon_start].unique()[0]
    except:
        pass
    try:
        max_c_val = final.c_val[final.chr_pos == exon_stop].unique()[0]
    except:
        pass
    final['nuc_change_tmp'] = final['nuc_change'].str.replace(':','>')
    #Add c values
    def conditions(final):
        if (final['AA_change'] == 'Intronic'):
            if (final['chr_pos'] < exon_start):   
                return min_c_val.__str__() + (final['chr_pos'] - exon_start).__str__()
            else:
                return max_c_val.__str__() +'+' + (final['chr_pos'] - exon_stop).__str__()
    final['c_val_tmp'] = final.apply(conditions, axis = 1)
    final['c_nom'] = 'NM_000059.4(BRCA2):c.'+ final['c_val_tmp'] + final['nuc_change_tmp']
    final['g_nom'] = 'NC_000013.11:g.' + final['chr_pos'].astype('str') + final['nuc_change_tmp']
    #Susbet splicing data
    splicing_exon = splicing[splicing.Exon == amplicon]
    splicing_exon = splicing_exon[['g.nom', 'Affect_splicing']]
    final = final.merge(splicing_exon, how = 'left', left_on = ['g_nom'], right_on = ['g.nom'])
    #Clean and save the file
    final.drop(['change', 'GRCh38Location', 'c_val', 'nuc_change_tmp', 'c_val_tmp'], axis = 1, inplace = True)
    final.to_excel('./variant_counts_no_thresh/annotated/'+ file_str +'_annotated.xlsx', na_rep = 'NA', index = False)
    print("Work done for " + i)
print("Work done for all files")




