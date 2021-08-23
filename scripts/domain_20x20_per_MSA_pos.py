
#read in (pileup) clustalo MSAs and make one 20x20 matrix (csv with subst counts) for each position in the alignment (across all sequences that are part of the MSA)

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
base_dir = '/storage1/hezscha/genome_proteome_map/results/domain_analysis/'
outdir = '/storage1/hezscha/genome_proteome_map/domain_analysis/results/'

#read in MSA: it has the file names and dom_start and dom_end
#for each seq, read in the prism_merged file and count up gnomad vars

#since we agreed w Amelie to just do one MSA lets get rid of the n
#df for counting subs
order_AA_WT = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'X']
order_AA_sub = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', '*']

#I'm using 0 (the first of ten alignments) for now, might change later
aln_file = os.path.join(base_dir, dom_name, 'aln_pileup', 'aln_0_pileup') 
len_aln = 0
#dict of dfs, one per position
mats_per_pos = {}
gaps_per_pos = {}
print(aln_file)

first = 1
for seq_record in SeqIO.parse(open(aln_file, mode='r'), 'fasta'):
	
	#all seqs in one alignment have the same length (since they are aligned to each other). When we read in the first seq, set up dfs for each position in mats_per_pos
	if first:
		for msa_pos,base in enumerate(seq_record.seq):
			mats_per_pos[msa_pos+1] = pd.DataFrame(index=order_AA_sub, columns=order_AA_WT).fillna(0)
			gaps_per_pos[msa_pos+1] = 0
			
		len_aln = len(seq_record.seq)
		first = 0
	
	#parse the record. in my case id = description = Name 
	ls = seq_record.name.split(';')
	uID = ls[0]
	dom_id = ls[1]
	dom_start = int(ls[2].split(':')[0])
	dom_end = int(ls[2].split(':')[1])
	if args.verbose:
		print(uID)
		print(dom_id, dom_start, dom_end)
		print(seq_record.seq)
	
	md, df = read_from_prism(os.path.join('/storage1/shared/data/prism_merge/', uID[0:2], uID[2:4], uID[4:6], 'prism_merged_'+uID+'.all.txt'))
	
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
	
	
	#msa_pos indexes the msa, prot_pos indexes the seq of the current protein.
	prot_pos = dom_start
	for msa_pos,base in enumerate(seq_record.seq):
		#add the gnomad vars at prot_pos to the sub mat at msa_pos
		#debug
		if msa_pos == 0:
			print(uID, base, msa_pos+1, prot_pos, list(df_gnom.loc[df_gnom['resi'] == prot_pos].variant))
		#debug end
		
		#if the current protein has a gap here, there can be no gnomad vars
		if base == '-':
			#can count gaps here if you want to
			gaps_per_pos[msa_pos+1] += 1
			#continue
			
		else:
			#debug: show the gnomad vars at this pos for checking
			# ~ if args.verbose:
				# ~ #if there are no gnomad vars at this pos
				# ~ if df_gnom.loc[df_gnom['resi'] == prot_pos].empty:
					# ~ print(base, msa_pos+1, prot_pos, 'no gnomad vars at this pos')
				# ~ else:	
					# ~ print(base, msa_pos+1, prot_pos, len(df_gnom.loc[df_gnom['resi'] == prot_pos]))
					# ~ #print(df_gnom.loc[df_gnom['resi'] == prot_pos])
			# ~ #debug end	
		
			#now add the gnomad vars at prot_pos to the sub mat at msa_pos
			#if there are no gnomad vars at this pos there is nothing to add
			if df_gnom.loc[df_gnom['resi'] == prot_pos].empty:
				pass
			else:	
				for index, row in df_gnom.loc[df_gnom['resi'] == prot_pos].iterrows():		
					mats_per_pos[msa_pos+1].loc[row['aa_var'][0],row['aa_ref'][0]] += 1	
			
			#only increase the position in the current protein if the character wasn't a gap
			#print('Increase prot_pos')
			prot_pos += 1

		
for msa_pos in range(1,len_aln+1):
	outfile = os.path.join(outdir, args.dom, "submat_"+str(msa_pos)+".csv")
	mats_per_pos[msa_pos].to_csv(outfile)	
