
#read in (pileup) clustalo MSAs and count the number of gnomad muts per position in the alignment (across all sequences that are part of the MSA)

import sys
import os
import requests
import argparse
import re
import numpy as np
import pandas as pd 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess

#local paths for HZ. Change this to your paths
#server	
path_base = '/storage1/hezscha/'
import_path_base = path_base + 'src/'

sys.path.insert(1, import_path_base + 'PRISM/software/domain_protein_features/scripts/')
from get_domains import extract_single_protein_pfam
sys.path.insert(1, import_path_base + 'helper_python_functions/')
from helpers import do_align,get_uniprot_seq,fasta_print
#sys.path.insert(1, import_path_base + 'gnomad_to_prism/scripts/')
#from vcf_parser_func_v09_revep import get_uniprot, get_prot_seq
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
#from PrismData_HZ_allow_SNP_v4 import PrismParser, VariantData
from PrismData import PrismParser, VariantData

def read_from_prism(primsfile):
	parser = PrismParser()
	dataframe = parser.read(primsfile).dataframe
	meta_data = parser.read_header(primsfile)
	return meta_data, dataframe

def write_prism(metadata, dataframe, prism_file, comment=''):
	variant_dataset = VariantData(metadata, dataframe)
	parser = PrismParser()
	parser.write(prism_file, variant_dataset, comment_lines=comment)

#parse arguments
################################################################################
#https://stackoverflow.com/questions/50021282/python-argparse-how-can-i-add-text-to-the-default-help-message
parser = argparse.ArgumentParser()
parser.add_argument('-dom', dest="dom", required = True, help='The domain to run for.')
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")
parser.add_argument('-n', dest="n",default=10, type = int, help="Number of permuted sequence lists to create.")
args = parser.parse_args()
################################################################################

dom_name = args.dom
base_dir = '/storage1/hezscha/genome_proteome_map/results/domain_analysis/'

#read in MSA: it has the file names and dom_start and dom_end
	#for each seq, read in the prism_merged file and count up gnomad vars

gnomad_vars_per_pos = {}
msa_lens = {}
for i in range(args.n):
	aln_file = os.path.join(base_dir, dom_name, 'aln_pileup', 'aln_'+str(i)+'_pileup') 
	len_aln = 0
	gnomad_vars_per_pos['msa'+str(i)] = {}
	print(aln_file)
	
	for seq_record in SeqIO.parse(open(aln_file, mode='r'), 'fasta'):
		#in my case id = description = Name 
		ls = seq_record.name.split(';')
		uID = ls[0]
		dom_id = ls[1]
		dom_start = int(ls[2].split(':')[0])
		dom_end = int(ls[2].split(':')[1])
		#print(uID)
		#print(dom_id, dom_start, dom_end)
		#print(seq_record.seq)
		
		if not len_aln:
			len_aln = len(seq_record.seq)
			if not len_aln in msa_lens:
				msa_lens[len_aln] = [i]
			else:
				msa_lens[len_aln].append(i)
		
		#we need the prism_merge file again
		md, df = read_from_prism(os.path.join('/storage1/shared/data/prism_merge/', uID[0:2], uID[2:4], uID[4:6], 'prism_merged_'+uID+'.all.txt'))
		#this should not fail because the corresponding prism_merge file was opened before to obtain the seq. Therefore no catch, if this fails I want to know
			
		#unlist the residue number
		df['resi'] = df['resi'].apply(lambda x: x[0])
		df['aa_ref'] = df['aa_ref'].apply(lambda x: x[0])
		df['aa_var'] = df['aa_var'].apply(lambda x: x[0])
		
		#again, we need to act once for each instance of the domain in the current protein
		#actually no, we should have one line per instance already in the alignment. We need to get the dom_id from the line
		#for dom_id in dom_dict[uID]:
		
		#dom_start = dom_dict[uID][dom_id]['dom_start_target']
		#dom_end = dom_dict[uID][dom_id]['dom_end_target']
		#msa_pos indexes the msa, prot_pos indexes the seq of the current protein.
		prot_pos = dom_start
		for msa_pos,base in enumerate(seq_record.seq):
			if not msa_pos+1 in gnomad_vars_per_pos['msa'+str(i)]:
				gnomad_vars_per_pos['msa'+str(i)][msa_pos+1] = 0
			
			#if the current protein has a gap here, there can be no gnomad vars
			if base == '-':
				continue
			else:
				#debug
				if args.verbose:
					if df.loc[df['resi'] == prot_pos].empty:
						print(prot_pos, 'empty df')
					else:	
						print(prot_pos, len(df.loc[df['resi'] == prot_pos]))
					#print(df.loc[df['resi'] == prot_pos])
				#debug
				#now count the number of vars at this position in current protein
				gnomad_vars_per_pos['msa'+str(i)][msa_pos+1] += len(df.loc[df['resi'] == prot_pos])
				#only increase the position inthe current protein if the character wasn't a gap
				prot_pos += 1

	gnom_file = os.path.join(base_dir, dom_name, 'gnom_counts', dom_name + '_gnomad_vars.csv_' + str(i)) 
	with open(gnom_file, 'w') as tabOUT:
		print('#MSA_pos nr_gomad_vars', file = tabOUT)
		for msa_pos in range(1,len_aln+1):
			print(msa_pos, gnomad_vars_per_pos['msa'+str(i)][msa_pos], file = tabOUT)

#now, average over all alignments with the same length
avg_gnom_vars = {}
for aln_len in msa_lens:
	#only avergage if there is more than one msa with this length
	if len(msa_lens[aln_len]) > 1:
		avg_gnom_vars[aln_len] = {}
		
		for msa_pos in range(1,aln_len+1):
			#make a list of all gnomad_var_count at this position across the msas.
			#I.e. if there are 4 MSAs with len 71, there will be 4 different counts for pos 1,2,3, .. 71
			all_vals = []
			#create a subdict for this pos
			avg_gnom_vars[aln_len][msa_pos] = {}
			for i in msa_lens[aln_len]:
				#use the saved counts: gnomad_vars_per_pos['msa'+str(i)]
				all_vals.append(gnomad_vars_per_pos['msa'+str(i)][msa_pos])
			
			#debug
			
			#debug
			avg_gnom_vars[aln_len][msa_pos]['avg'] = np.mean(all_vals)
			avg_gnom_vars[aln_len][msa_pos]['std'] = np.std(all_vals)
			avg_gnom_vars[aln_len][msa_pos]['median'] = np.median(all_vals)
	
#go through and print the averaged gnomad var counts as well		
for aln_len in avg_gnom_vars:
	print(aln_len)
	gnom_file = os.path.join(base_dir, dom_name, 'gnom_counts', dom_name + '_gnomad_vars.csv_avg_' + str(aln_len))
	with open(gnom_file, 'w') as tabOUT:
		print('#len:', aln_len, 'nr MSAs:', len(msa_lens[aln_len]), file = tabOUT)
		print('#MSA_pos avg_gomad_vars std_gomad_vars median_gomad_vars', file = tabOUT)
		for msa_pos in range(1,aln_len+1):
			print(msa_pos, avg_gnom_vars[aln_len][msa_pos]['avg'], avg_gnom_vars[aln_len][msa_pos]['std'], avg_gnom_vars[aln_len][msa_pos]['median'], file = tabOUT)
	
