
#only works on nike with
#conda activate p3-env-cyvcf2
#because it invokes clustalo as a subprocess which is installed into that env

#read in a fasta file with protein sequence for one domain (created by another script) and create n permuted/shuffled versions of the list, then align them with clustalo

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
import random

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

print('Starting shuffle.\nEnsure you are on nike and have activated conda env p3-env-cyvcf2.')

dom_name = args.dom
base_dir = '/storage1/hezscha/genome_proteome_map/results/domain_analysis/'

if not os.path.exists(os.path.join(base_dir, dom_name)):
	print('There is no directory for', dom, '. Perhaps need to run make_domain_fastas.py first.')
	sys.exit()

if not os.path.exists(os.path.join(base_dir, dom_name, 'shuffle')): 
	os.makedirs(os.path.join(base_dir, dom_name, 'shuffle'))
if not os.path.exists(os.path.join(base_dir, dom_name, 'aln')): 
	os.makedirs(os.path.join(base_dir, dom_name, 'aln'))
	
seq_file = os.path.join(base_dir, dom_name, dom_name + '_seqs.fasta')
fasta_list = []
for seq_record in SeqIO.parse(open(seq_file, mode='r'), 'fasta'):
	fasta_list.append(seq_record)

for i in range(args.n):
	print('iteration:', i)
	f_out = os.path.join(base_dir, dom_name, 'shuffle', dom_name + '_seqs.fasta_' + str(i))
	aln_file = os.path.join(base_dir, dom_name, 'aln', dom_name + '_aligend.fasta_' + str(i)) 
	print('Shuffling list...')
	b = random.sample(fasta_list, len(fasta_list))
	
	#print shuffled seq file
	with open(f_out, "w") as output_handle:
		r=SeqIO.write(b, output_handle, 'fasta')
	
	print('Making alignment...')
	#align shuffled seq file. Will be slow because I'm using pileup
	shell_command = ['clustalo', '-i', f_out, '-o', aln_file, '--force', '--pileup']
	subprocess.run(shell_command,stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)



