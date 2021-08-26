
#bascially for Amelie's DHFR question
#similar to domain_count_gnomad_vars_per_MSA_pos_v3.py, but also counting clinvar_patho and benign and omitting skyling for now (can be added in later).
#read in (pileup) clustalo MSAs and count the number of gnomad muts per position in the alignment (across all sequences that are part of the MSA)
#also report the number of gaps per MSA position and the information content calculated by skylign

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
args = parser.parse_args()
################################################################################

dom_name = args.dom
base_dir = '/storage1/hezscha/genome_proteome_map/domain_analysis/results/'


scores_per_pos = {}

#gnomad vars and gaps:
#read in MSA: it has the file names and dom_start and dom_end
	#for each seq, read in the prism_merged file and count up gnomad vars

aln_file = os.path.join(base_dir, dom_name, 'alignment', 'aln_seed_instance') 
len_aln = 0
print(aln_file)

for seq_record in SeqIO.parse(open(aln_file, mode='r'), 'fasta'):
	
	#all seqs in one alignment have the same length (since they are aligned to each other). When we read in the first seq, set up dfs for each position in mats_per_pos
	if first:
		for msa_pos,base in enumerate(seq_record.seq):
			scores_per_pos[msa_pos+1] = {}
			scores_per_pos[msa_pos+1]['gnom_vars'] = 0
			scores_per_pos[msa_pos+1]['n_gaps'] = 0
			scores_per_pos[msa_pos+1]['bits'] = 0
			
			# ~ gaps_per_pos[msa_pos+1] = 0
			# ~ WT_pfam_per_pos[msa_pos+1] = {}
			# ~ WT_human_per_pos[msa_pos+1] = {}
			# ~ for base in order_AA_WT+['-']:
				# ~ WT_pfam_per_pos[msa_pos+1][base] = 0
				# ~ WT_human_per_pos[msa_pos+1][base] = 0	
			
		len_aln = len(seq_record.seq)
		first = 0
	
	#parse the record. in my case id = description = Name 
	#pfam seq: not counting these
	if '/' in seq_record.name:
			continue
	
	#human proteome seq, count gnomad vars 
	else:
		
		#in my case id = description = Name 
		ls = seq_record.name.split(';')
		uID = ls[0]
		dom_id = ls[1]
		dom_start = int(ls[2].split(':')[0])
		dom_end = int(ls[2].split(':')[1])
		if args.v > 0:
			print(uID)
			print(dom_id, dom_start, dom_end)
			print(seq_record.seq)
		
		#we need the prism_merge file again
		md, df = read_from_prism(os.path.join('/storage1/shared/data/prism_merge/', uID[0:2], uID[2:4], uID[4:6], 'prism_merged_'+uID+'.all.txt'))
		#this should not fail because the corresponding prism_merge file was opened before to obtain the seq. Therefore no catch, if this fails I want to know
			
		#if there is no gnomad component file, skip to next
		#for else construct: if we didn't break the inner for loop, skip to next iteration of the outer for loop with continue
		for File in md['merged']:
			if 'prism_gnomad_001' in md['merged'][File]:
				gnom_file_nr = '_'+File.split('_')[-1]
				break
		else:
			print(seq_record.name, 'has no gnomad component in the merge file')
			continue	
			
		#unlist the residue number
		df['resi'] = df['resi'].apply(lambda x: x[0])
		df['aa_ref'] = df['aa_ref'].apply(lambda x: x[0])
		df['aa_var'] = df['aa_var'].apply(lambda x: x[0])
		
		#subset the df to only lines that have a gnomad AF
		df_gnom = df.loc[df['gnomad;AF_tot'+gnom_file_nr].notnull()]
		
		#again, we need to act once for each instance of the domain in the current protein
		#actually no, we should have one line per instance already in the alignment. We need to get the dom_id from the line
		#for dom_id in dom_dict[uID]:
		
		#dom_start = dom_dict[uID][dom_id]['dom_start_target']
		#dom_end = dom_dict[uID][dom_id]['dom_end_target']
		#msa_pos indexes the msa, prot_pos indexes the seq of the current protein.
		prot_pos = dom_start
		for msa_pos,base in enumerate(seq_record.seq):
			#if not msa_pos+1 in scores_per_pos:
			#	scores_per_pos[msa_pos+1] = 0
			
			#if the current protein has a gap here, there can be no gnomad vars
			if base == '-':
				scores_per_pos[msa_pos+1]['n_gaps'] += 1
				continue
			else:
				#debug
				if args.verbose:
					print(uID, base, msa_pos+1, prot_pos, list(df_gnom.loc[df_gnom['resi'] == prot_pos].variant))
					# ~ if df.loc[df['resi'] == prot_pos].empty:
						# ~ print(prot_pos, 'empty df')
					# ~ else:	
						# ~ print(prot_pos, len(df.loc[df['resi'] == prot_pos]))
					#print(df.loc[df['resi'] == prot_pos])
				#debug
				
				#now count the number of vars at this position in current protein
				scores_per_pos[msa_pos+1]['gnom_vars'] += len(df_gnom.loc[df['resi'] == prot_pos])
				#only increase the position inthe current protein if the character wasn't a gap
				prot_pos += 1

#read in bits file from skylign to get the information content of the position, i.e. column, which is the sum of the information content of all 20 aas in that column.

outfile_file = os.path.join(base_dir, dom_name, 'gnom_counts', 'scores.csv') 
with open(outfile_file, 'w') as tabOUT:
	print('#MSA_pos nr_gomad_vars nr_gaps bits', file = tabOUT)
	for msa_pos in range(1,len_aln+1):
		print(msa_pos, scores_per_pos[msa_pos]['gnom_vars'], scores_per_pos[msa_pos+1]['n_gaps'], scores_per_pos[msa_pos+1]['bits'], file = tabOUT)


