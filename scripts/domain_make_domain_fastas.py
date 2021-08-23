
#go through a list of proteins and extract the asked for domain into a fasta
#the protein list should already be made by looking for prots with that domain, i.e. with stats_make_domain_lists.py

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
parser.add_argument('-uniprot', dest="uniprot", help="The uniprot ID for which to merge files. The authorative sequence for merging is selected with following hierarchy: 1. transcript with most clinvar vars 2. transcript with most gnomad vars 3. transcript with most spliceAI vars 4. uniprot isoform 1, unless otherwise specified with -seq")
parser.add_argument('-dom', dest="dom", required = True, help='The domain to run for.')
parser.add_argument('-ff','--fromfile', dest="ff", help="A file from which to read uniprot IDs (one per line).")
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")
args = parser.parse_args()
################################################################################

dom_name = args.dom
base_dir = '/storage1/hezscha/genome_proteome_map/results/domain_analysis/'

#go through all merge files
#get domain info directly from pfam 
#pfam_dict = extract_single_protein_pfam( uniprot_id, verbose=False )
#Somehow map the domain info (which is relative to uniprot isoform 1) to the target seq (the one in the metadata). Probably via alignment of the uniprot 1 and target seqs. Check if they are identical before

if args.uniprot:
	uniprot_id_list = args.uniprot.split(',')
elif args.ff:
	uniprot_id_list = []
	with open(args.ff, 'r') as IN:
		for line in IN:
			uniprot_id_list.append(line.rstrip())

#record the start and end of domains in the target_seq
dom_dict = {}

if not os.path.exists(os.path.join(base_dir, dom_name)):
	os.makedirs(os.path.join(base_dir, dom_name))
	
OUT = open(os.path.join(base_dir, dom_name, dom_name + '_seqs.fasta'), 'w')
for uID in uniprot_id_list:
	print(uID)
	try:
		md, df = read_from_prism(os.path.join('/storage1/shared/data/prism_merge/', uID[0:2], uID[2:4], uID[4:6], 'prism_merged_'+uID+'.all.txt'))
	except:
		print("Couldn't open prism_merge file for", uID)
		continue
	
	dom_dict[uID] = {}
	#uID = 'A1IGU5'
	target_seq = md['protein']['sequence']
	uniprot_seq = get_uniprot_seq(uID)

	#get domain mapping to uniprot isoform 1
	pfam_dict = extract_single_protein_pfam(uID, verbose=False)

	#transfer domain mapping to target seq
	if not uniprot_seq == target_seq:
		align_obj = do_align(uniprot_seq, target_seq)
		
		#hack for now:
		i = 1
		for dom in pfam_dict:
			if dom['id'] == dom_name:
				dom_id = dom_name+'_'+str(i)
				i += 1
				
				dom_dict[uID][dom_id] = {}
				dom_dict[uID][dom_id]['dom_start_uni'] = dom['start']
				dom_dict[uID][dom_id]['dom_end_uni'] = dom['end']
				if args.verbose:
					print(dom_dict[uID][dom_id]['dom_start_uni'], dom_dict[uID][dom_id]['dom_end_uni'])
				#break
		
		if args.verbose > 1:
			print(dom_dict[uID])
		#debug: skip rest for now so I can see if the dom_dict is correct
		#continue
		
		
		#we need to identify which pos in target_seq matches to the start and end of the domain
		#see example_alignment_case
		#we need dom_start_target and dom_end_target, the start and stop of the domain in the (unaligned) target seq
		#I see that:
		#dom_start_uni_aln = dom_start_uni + n_gaps_uni_before
		#and dom_start_target_aln == dom_start_uni_aln
		#and dom_start_target = dom_start_target_aln - n_gaps_target_before
		#need n_gaps_uni_before_start, n_gaps_uni_before_end, n_gaps_target_before_start and n_gaps_target_before_end
		
		uniprot_seq_aln = align_obj.aligned_a
		target_seq_aln = align_obj.aligned_b
		
		#we need to act on every instance of the domain found in this protein.
		for dom_id in dom_dict[uID]:
			print(dom_id)
		
			#first, we count all gaps up to dom_start_uni:
			n_gaps_uni_before_start = uniprot_seq_aln[:dom_dict[uID][dom_id]['dom_start_uni']-1].count('-')
			#then, we check if the next character is a non-gap. If not, n_gaps_uni_before++
			for pos,l in enumerate(uniprot_seq_aln[dom_dict[uID][dom_id]['dom_start_uni']-1:]):
				if l != '-':
					break
				else:
					n_gaps_uni_before_start += 1
					
			#gaps before the end of domain:
			n_gaps_uni_before_end = uniprot_seq_aln[:dom_dict[uID][dom_id]['dom_end_uni']-1].count('-')
			#then, we check if the next character is a non-gap. If not, n_gaps_uni_before++
			for pos,l in enumerate(uniprot_seq_aln[dom_dict[uID][dom_id]['dom_end_uni']-1:]):
				if l != '-':
					break
				else:
					n_gaps_uni_before_end += 1		
			
			#now we can calc where the start of the domain is in the aligned uniprot seq. This is also the start of the domain in the aligned target_seq. (the premise of alignment is that pos with the same number are equivalent between the two seqs) 
			dom_start_uni_aln = dom_dict[uID][dom_id]['dom_start_uni'] + n_gaps_uni_before_start
			dom_end_uni_aln = dom_dict[uID][dom_id]['dom_end_uni'] + n_gaps_uni_before_end
			dom_start_target_aln = dom_start_uni_aln
			dom_end_target_aln = dom_end_uni_aln
				
			#now that we have dom_start_target/uni_aln, we just count gaps up to there	
			#gaps before the domain in target_seq
			n_gaps_target_before_start = target_seq_aln[:dom_start_target_aln-1].count('-')		
			#gaps before the end of domain:
			n_gaps_target_before_end = target_seq_aln[:dom_end_target_aln-1].count('-')
						
			#so we again calculate the start and end of the domain in the unaligned target_seq, now that we have the gap counts:
			dom_dict[uID][dom_id]['dom_start_target'] = dom_start_target_aln - n_gaps_target_before_start
			dom_dict[uID][dom_id]['dom_end_target'] = dom_end_target_aln - n_gaps_target_before_end
	 
			#debug
			if args.verbose:
				print(dom_dict[uID][dom_id]['dom_start_target'], dom_dict[uID][dom_id]['dom_end_target'])
			#debug
	 
			#print the chunk of the aligned target_seq at the place where the domain is (may inlcude gaps)
			print('>', uID, ';', dom_id, ';', dom_dict[uID][dom_id]['dom_start_target'], ':', dom_dict[uID][dom_id]['dom_end_target'], sep = '', file = OUT)
			fasta_print(target_seq[dom_dict[uID][dom_id]['dom_start_target']-1:dom_dict[uID][dom_id]['dom_end_target']], OUT)
		
	#else, if uniprot_seq == target_seq we can directly use the positions given by pfam 
	else:
		if args.verbose:
			print('Same')
		#hack for now:
		i = 1
		for dom in pfam_dict:
			if dom['id'] == dom_name:
				dom_id = dom_name+'_'+str(i)
				i += 1
				
				dom_dict[uID][dom_id] = {}
				dom_dict[uID][dom_id]['dom_start_target'] = dom['start']
				dom_dict[uID][dom_id]['dom_end_target'] = dom['end']
				if args.verbose:
					print(dom_dict[uID][dom_id]['dom_start_target'], dom_dict[uID][dom_id]['dom_end_target'])
				#break
		
		if args.verbose > 1:		
			print(dom_dict[uID])
		#debug: skip rest for now so I can see if the dom_dict is correct
		#continue		
		
		#we need to act on every instance of the domain found in this protein.
		for dom_id in dom_dict[uID]:
			print(dom_id)
		
			print('>', uID, ';', dom_id, ';', dom_dict[uID][dom_id]['dom_start_target'], ':', dom_dict[uID][dom_id]['dom_end_target'], sep = '', file = OUT)
			fasta_print(target_seq[dom_dict[uID][dom_id]['dom_start_target']-1:dom_dict[uID][dom_id]['dom_end_target']], OUT)

OUT.close()
