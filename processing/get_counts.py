#!/usr/bin/env python
#coding: utf-8
#Script to get count of reads in each sam files
#To be run in directory with sam files
#Load libraries
import sys
import os
import subprocess
import operator
import pickle
from Bio.Seq import Seq
import numpy as np
import re
import string
#Get list of amplicons by checking the sam files available
all_files = os.listdir()
sam_list = [x for x in all_files if '.sam' in x]
amplicon_list = [s.split('_')[0] for s in sam_list]
amplicon_list = list(np.unique(amplicon_list))

#Group experiments and file names
exp_groupings_str = ''
for amp in amplicon_list:
    for replicate in ['RP1', 'RP2']:
        for comp in ['OLA', 'CIS', 'DMSO']:
            exp_groupings_str = exp_groupings_str + amp + replicate + comp + ','
            exp_groupings_str = exp_groupings_str + amp + '_PRE' + replicate + ','
            exp_groupings_str = exp_groupings_str + amp + '_' + comp + replicate + '+'
exp_groupings_str = exp_groupings_str[:-1]
exp_groupings_list = exp_groupings_str.split('+')
exp_groupings = {}

for grouping in exp_groupings_list:
    g_info = grouping.split(',')
    exp_groupings[g_info[0]] = g_info[1:]
    
#Create a dictionary with reference for each exon
my_dir =  os.getcwd()
fasta_dir = '../../../../../reference'
ref_seqs = {}
for amplicon in amplicon_list:
    ref_file = open(fasta_dir + '/' + amplicon + '.fa', 'r')
    ref_header = ref_file.readline()
    ref_seq = ref_file.readline().strip()
    ref_seqs[amplicon] = ref_seq
    ref_file.close()
    
#Load edits info for each amplicon 
edits_per_amp = {}
editing_info = open('../../../../../reference/BRCA2_editing_data.txt', 'r')
editing_info_header = editing_info.readline()
for line in editing_info:
	edits = line.strip().split()
	edits_per_amp[edits[0]] = edits[1:]
    
#Sets a reads threshold for making it into the output
alignment_score_threshold = 200

#Function to convert cigar to edit string
def cigar_to_edits(cigar, seq, ref_seq):
	MDI_indexes = []
	for x in range(0,len(cigar)):
		if (cigar[x] == 'M') or (cigar[x] == 'D') or (cigar[x] == 'I'):
			MDI_indexes.append(x)
	cigar_index = 0
	read_index = 0
	ref_index = 0
	mod_cigar = ''
	edit_string = ''
	for mdi_index in MDI_indexes:
		char = cigar[mdi_index]
		length = int(cigar[cigar_index:mdi_index])
		cigar_index = mdi_index+1
		if char == 'M':
			match_length = 0
			adj = ref_index-read_index
			for x in range(ref_index,ref_index+length):
				if seq[x-adj] == ref_seq[x]:
					match_length+=1
					if x == (ref_index+length):
						mod_cigar += (str(match_length)+'M')
						match_length = 0
				elif seq[x-adj] != ref_seq[x]:
					if match_length >= 1:
						mod_cigar += (str(match_length)+'M'+'1X'+seq[x-adj])
					elif match_length == 0:
						mod_cigar += ('1X'+seq[x-adj])
					edit_string += (str(ref_index+x+1)+'-'+'X-'+seq[x-adj]+',')
					match_length = 0
			#adjust the read_index and ref_index
			ref_index += length
			read_index += length
		elif char == 'D':
			mod_cigar += (str(length)+'D')
			edit_string += (str(ref_index+1)+'-D'+str(length)+',')
			ref_index += length
		elif char == 'I':
			insertion_sequence = seq[read_index:read_index+length]
			mod_cigar += (str(length)+'I'+insertion_sequence)
			edit_string += (str(ref_index+1)+'-I'+str(length)+'-'+insertion_sequence+',')
			read_index += length
	if edit_string == '':
		edit_string = '-WT'
	else:
		edit_string = edit_string[:-1]
	return([mod_cigar, edit_string])

#dol_lookup -- a function to lookup a key-number pairing in a dicitonary -- instead of returning a key error, return a 0 if not present)
def dol_lookup(counts_dol,key):
	if key in counts_dol:
		return counts_dol[key][0]
	else:
		return 0

#returns a dictionary for all variants
def sam_to_variant_cigar_dict(sam_file, min_alignment_score):
	line_count = -2
	variant_dict = {}
	reads_not_aligning = 0
	for line in sam_file:
		if line_count < 0:
			line_count += 1
			continue
		else:
			line_count += 1
			sam_data = line.strip().split('\t')
			cigar = sam_data[5]
			AS = sam_data[11]
			seq = sam_data[9]
			index_of_score = AS.find('i:')+2
			score = float(AS[index_of_score:])
			if score > min_alignment_score:
				if seq in variant_dict:
					variant_dict[seq][0]+=1
				else:
					variant_dict[seq]=[1,cigar]
			else: 
				reads_not_aligning += 1
	print(reads_not_aligning, "reads not aligning in sam file ", str(sam_file), " out of ", line_count, " total reads.") 
	return variant_dict, line_count
    
#Sort count dictionary
def sort_counts_dol(counts_dol):
    return sorted(counts_dol.keys(), key=lambda k: -1*int(counts_dol[k][0]))

#Return total read counts of dictionary
def get_read_count_dol(counts_dol):
    total_reads = 0
    for my_key in counts_dol:
        total_reads += counts_dol[my_key][0]
    return total_reads

# Create all possible SNV for a sequence
def mutagenize(my_string):
	upper_string = my_string.upper()
	all_variants = [upper_string]
	for i in range (0,len(upper_string)):
		for j in ('A', 'C', 'G', 'T'):
			if j != upper_string[i]:
				all_variants += [upper_string[:i]+j+upper_string[i+1:]]
	return all_variants

# Get the list of all HDR variants + WT + PAM
def get_all_hdr_variants(amplicon,ref_seq_dict,edits_per_amp_dict): #requires mutagenize function above
    ref_seq = ref_seq_dict[amplicon]
    all_hdr_variants = [ref_seq] #list will start with total wt
    my_editing_info = edits_per_amp_dict[amplicon]
    edits5 = my_editing_info[0].strip().split(',')
    edits3 = my_editing_info[1].strip().split(',')
    mut_start = int(my_editing_info[4])
    mut_end = int(my_editing_info[5])
    all_pam_edits = edits5 + edits3
    mut_ref_seq = str(ref_seq)
    for edit in all_pam_edits:
        edit_location = int(edit[:edit.find('-')])
        new_base = edit[-1]
        mut_ref_seq = mut_ref_seq[:edit_location-1]+new_base+mut_ref_seq[edit_location:]
    f_adapt = mut_ref_seq[:mut_start-1]
    mut_seq = mut_ref_seq[mut_start-1:mut_end]
    r_adapt = mut_ref_seq[mut_end:]
    all_mut_seqs = mutagenize(mut_seq)
    for seq in all_mut_seqs:
        all_hdr_variants.append(f_adapt+seq+r_adapt)
    return all_hdr_variants

#Dictionary with codons/AA
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
    'TAA': '_', 'TAG': '_', 'TGA': '_'
}

#Add synonymous/non-synonymous/non-sense information
def seq_to_codons(ref_seq, edit, edits3, edits5, seq, genomic_pos, exon_start):
    if edit == '-WT' or edit == edits3[0] or edit == edits5[0] or edit == edits3[0]+','+edits5[0] or edit == edits5[0]+','+edits3[0]:
        coding = 'Synonymous'
    else: 
        ref_seq = list(ref_seq)
        edit = re.sub(edits3[0], '', edit)
        edit = re.sub(edits5[0], '', edit)
        edit = re.sub(',', '', edit)
        pos = int(edit.split('-')[0])
        ref_seq_var = list(seq)
        ref_seq_var[pos - 1] = ref_seq[pos - 1]
        frame = (exon_start - genomic_pos) % 3 
        ref_seq_var = ''.join(ref_seq_var)
        ref_seq_var = ref_seq_var[frame:]
        seq = seq[frame:]
        coding = 'Synonymous'
        for i in range(0, len(seq),3):
            codon_seq = seq[i:i+3]
            if len(codon_seq) < 3 or len(ref_seq_var[i:i+3]) < 3 :
                break       
            codon = DNA_Codons[codon_seq]
            seq_codon_ref_seq = ref_seq_var[i:i+3]
            codon_ref_seq = DNA_Codons[seq_codon_ref_seq]
            if codon == codon_ref_seq:
                pass
            elif codon == '_':
                coding = 'Non-sense'
            else:
                coding = 'Non-synonymous'
        return coding
    
#Get the sequence for the variant with only the 3' PAM   
def get_3_edit_only_variant(ref_seq, variant, edit, edits3, edits5):
    edit = edit.split(',')
    edit.remove(edits3[0])
    edit.remove(edits5[0])
    edit = edit[0]
    pos = int(edit.split('-')[0])
    nuc_changed = edit.split('-')[2]
    pos_edit_3 = int(edits3[0].split('-')[0])
    nuc_changed_edit_3 = edits3[0].split('-')[2]
    ref_seq_list = list(ref_seq)
    ref_seq_list[pos-1] = nuc_changed
    ref_seq_list[pos_edit_3-1] = nuc_changed_edit_3
    var_3_edited = ''.join(ref_seq_list)
    return(var_3_edited)

#Get the sequence for the variant with only the 5' PAM
def get_5_edit_only_variant(ref_seq, variant, edit, edits3, edits5):
    edit = edit.split(',')
    edit.remove(edits3[0])
    edit.remove(edits5[0])
    edit = edit[0]
    pos = int(edit.split('-')[0])
    nuc_changed = edit.split('-')[2]
    pos_edit_5 = int(edits5[0].split('-')[0])
    nuc_changed_edit_5 = edits5[0].split('-')[2]
    ref_seq_list = list(ref_seq)
    ref_seq_list[pos-1] = nuc_changed
    ref_seq_list[pos_edit_5-1] = nuc_changed_edit_5
    var_5_edited = ''.join(ref_seq_list)
    return(var_5_edited)

#Get each alternative version of the PAM sequence (only 3, only 5 and reverse)
def get_3_5_edit_PAM(ref_seq, edits3, edits5):
    pos3 = int(edits3[0].split('-')[0])
    pos5 = int(edits5[0].split('-')[0])
    nuc_changed_edit_3 = edits3[0].split('-')[2]
    nuc_changed_edit_5 = edits5[0].split('-')[2]
    PAM3 = list(ref_seq)
    PAM5 = list(ref_seq)
    PAM3[pos3-1] = nuc_changed_edit_3
    PAM5[pos5-1] = nuc_changed_edit_5
    PAM3 = ''.join(PAM3)
    PAM5 = ''.join(PAM5)
    seqPAM3 = Seq(PAM3)
    seqPAM5 = Seq(PAM5)
    revPAM3 = seqPAM3.reverse_complement()
    revPAM5 = seqPAM5.reverse_complement()
    return(PAM3, PAM5, revPAM3, revPAM5)

#Create the dictionary with all variants
working_dir = os.getcwd()
dict_of_var_dicts = {}
sample_list = []
total_reads_dict = {}

#Counts all reads in each sam file
print("Starting working on the sam files. This can take some time")
for i in os.listdir(os.getcwd()):
	if not i.endswith('.sam'):
		continue
	elif i.endswith('.sam'):
		index_of_first_dot = i.find('.')
		sample = i[:index_of_first_dot]
		index_of_first_under = i.find('_')
		timepoint = i[index_of_first_under+1:index_of_first_dot]
		experiment = i[:index_of_first_under]
		if sample not in sample_list:
			sample_list.append(sample)
		with open(i, 'r') as sam_file:
			dict_of_var_dicts[sample], total_reads_dict[sample] = sam_to_variant_cigar_dict(sam_file,alignment_score_threshold)            

#Parse the variants dictionary to get each count
header_string = 'variant\tcigar\text_cigar\tedit_string\tpre\tpost\tpre_freq\tpost_freq\tcodon'
for exp in exp_groupings.keys():
    if 'R' in exp:
        amplicon = exp[:exp.find('R')]
        print('Start working on ', amplicon)
    else:
        amplicon = exp
    ref_seq = ref_seqs[amplicon]
    my_editing_info = edits_per_amp[amplicon]
    total_reads = get_read_count_dol(dict_of_var_dicts[exp_groupings[exp][1]])
    edits5 = my_editing_info[0].strip().split(',')
    edits3 = my_editing_info[1].strip().split(',')
    genomic_pos = int(my_editing_info[6])
    exon_start = int(my_editing_info[8])
    exon_stop = int(my_editing_info[9])
    total_pre = total_reads_dict[exp_groupings[exp][0]]
    total_post = total_reads_dict[exp_groupings[exp][1]]
    #Get the dictionary with all possibles HDR variants
    all_hdr_variants = get_all_hdr_variants(amplicon, ref_seqs, edits_per_amp)
    #Extract counts for all variants in the pre samples
    pre_var_dict = dict_of_var_dicts[exp_groupings[exp][0]]
    post_var_dict = dict_of_var_dicts[exp_groupings[exp][1]]
    pre_var_rc = float(get_read_count_dol(pre_var_dict))
    post_var_rc = float(get_read_count_dol(post_var_dict))
    #Create a list of all variants present in either pre_var_dict or post_var_dict
    pre_variants = sort_counts_dol(pre_var_dict)
    post_variants = sort_counts_dol(post_var_dict)
    unique_post_variants = []
    for post_variant in post_variants:
        if post_variant not in pre_var_dict:
            unique_post_variants.append(post_variant)
    pre_post_variants = pre_variants+unique_post_variants
    exp_variants = pre_post_variants
    exp_master_dict = {}
    #Try to find the cigar for WT. If not found, we just need to look for the reverse
    try:
        cigar = pre_var_dict[all_hdr_variants[0]][1]
    except:
        seq = Seq(all_hdr_variants[0])
        reverse = seq.reverse_complement()
        cigar = pre_var_dict[str(reverse)][1]
    #Start working on each variant in the list of all possible HDR variants
    for variant in all_hdr_variants:
        seq = Seq(variant)
        reverse = seq.reverse_complement()
        outline_data = []
        edit = cigar_to_edits(cigar,variant,ref_seq)[1]
        outline_data += [cigar]
        outline_data += cigar_to_edits(cigar,variant,ref_seq)
        #Add count for the sequence + the reverse
        raw_reads = [dol_lookup(pre_var_dict,variant),dol_lookup(post_var_dict,variant)]
        raw_reads[0] += dol_lookup(pre_var_dict,reverse)
        raw_reads[1] += dol_lookup(post_var_dict,reverse)
        #For experiment with two PAM, add info counts for variants with only one PAM
        if edits3 != edits5 and edit != '-WT' and edit != edits3[0] and edit != edits5[0] and edit != edits3[0]+','+edits5[0] and edit != edits5[0]+','+edits3[0]:
            try:
                var_3_edited = get_3_edit_only_variant(ref_seq,variant, edit, edits3, edits5)
            except:
                print(edit, "NA")
            seq3 = Seq(var_3_edited)
            reverse3 = seq3.reverse_complement()          
            raw_reads[0] += dol_lookup(pre_var_dict,var_3_edited)
            raw_reads[0] += dol_lookup(pre_var_dict,reverse3)
            raw_reads[1] += dol_lookup(post_var_dict,var_3_edited)
            raw_reads[1] += dol_lookup(post_var_dict,reverse3)
            try:
                var_5_edited = get_5_edit_only_variant(ref_seq,variant, edit, edits3, edits5)
            except:
                print(edit, "NA")
            seq5 = Seq(var_5_edited)
            reverse5 = seq5.reverse_complement()  
            raw_reads[0] += dol_lookup(pre_var_dict,var_5_edited)
            raw_reads[0] += dol_lookup(pre_var_dict,reverse5)
            raw_reads[1] += dol_lookup(post_var_dict,var_5_edited)
            raw_reads[1] += dol_lookup(post_var_dict,reverse5)
        #Save for values just the PAM
        if edits3 != edits5 and (edit == edits3[0] + ',' + edits5[0] or edit == edits5[0] + ',' + edits3[0]):
            PAM3, PAM5, revPAM3, revPAM5 = get_3_5_edit_PAM(ref_seq, edits3, edits5)
            raw_reads[0] += dol_lookup(pre_var_dict,PAM3)
            raw_reads[0] += dol_lookup(pre_var_dict,PAM5)
            raw_reads[0] += dol_lookup(pre_var_dict,revPAM3)
            raw_reads[0] += dol_lookup(pre_var_dict,revPAM5)
            raw_reads[1] += dol_lookup(post_var_dict,PAM3)
            raw_reads[1] += dol_lookup(post_var_dict,PAM5)
            raw_reads[1] += dol_lookup(post_var_dict,revPAM3)
            raw_reads[1] += dol_lookup(post_var_dict,revPAM5)
        #Compute frequencies
        pre_freq = raw_reads[0] / total_pre
        post_freq = raw_reads[1] / total_post
        #Add data
        outline_data += raw_reads 
        outline_data.append(pre_freq)
        outline_data.append(post_freq) 
        change_in_codon = seq_to_codons(ref_seq, edit, edits3, edits5, variant, genomic_pos, exon_start)
        outline_data += [change_in_codon]
        exp_master_dict[variant] = outline_data
    #Save
    outfile_path = my_dir +'/variant_counts_no_thresh/' + exp + '.txt'
    outfile_open = open(outfile_path,'w')
    outfile_open.write(header_string+'\n')
    for variant in all_hdr_variants:
        output_list = [variant]+[str(i) for i in exp_master_dict[variant]]
        outfile_open.write('\t'.join(output_list) + '\n')
    print('Finished creating output file for '+ exp +' experiment.')
    outfile_open.close()



