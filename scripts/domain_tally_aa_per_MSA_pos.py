#read in the MSA of the pileup (all human domain instances) and the pfam seed alignment. for each MSA pos, tally how many of each WT aa occur at each pos, separately for seqs from human proteome (all human domain instances) and seqs from pfam seed

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

order_AA_WT = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'X', '-']

aln_file = os.path.join(base_dir, dom_name, 'pfam_to_instances', 'aln_seed_instance') 
len_aln = 0
#dict of dicts, one per position
WT_human_per_pos = {}
WT_pfam_per_pos = {}
print(aln_file)

first = 1
for seq_record in SeqIO.parse(open(aln_file, mode='r'), 'fasta'):
	
	#all seqs in one alignment have the same length (since they are aligned to each other). When we read in the first seq, set up dfs for each position in mats_per_pos
	if first:
		for msa_pos,base in enumerate(seq_record.seq):
			WT_pfam_per_pos[msa_pos+1] = {}
			WT_human_per_pos[msa_pos+1] = {}
			for base in order_AA_WT:
				WT_pfam_per_pos[msa_pos+1][base] = 0
				WT_human_per_pos[msa_pos+1][base] = 0
			
		len_aln = len(seq_record.seq)
		first = 0
	
	#parse the record. in my case id = description = Name 
	#pfam seq
	if '/' in seq_record.name:
		for msa_pos,base in enumerate(seq_record.seq):
			WT_pfam_per_pos[msa_pos+1][base] += 1
	
	#human proteome seq
	else:
		for msa_pos,base in enumerate(seq_record.seq):
			WT_human_per_pos[msa_pos+1][base] += 1
		

#once for pfam seed alignment
for msa_pos in range(1,len_aln+1):
	outfile = os.path.join(outdir, args.dom, "AA_distro_human_proteome_"+str(msa_pos)+".txt")
	with open(outfile, 'w') as OUT:
		for base in order_AA_WT:
			print(base, WT_human_per_pos[msa_pos][base], file = OUT)
				
#once for human proteome instances	
for msa_pos in range(1,len_aln+1):
	outfile = os.path.join(outdir, args.dom, "AA_distro_pfam_seed_"+str(msa_pos)+".txt")
	with open(outfile, 'w') as OUT:
		for base in order_AA_WT:
			print(base, WT_pfam_per_pos[msa_pos][base], file = OUT)
		 
	
