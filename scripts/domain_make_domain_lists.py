
#go through all prism_uniprot files and make lists of which domains occur in which protein. The lists are per domain and list all proteins this domain occurs in. The domains are named with their pfam name
#Just the domains occuring in the most proteins. So if domain A occurs in p1, p2 and p3, count 3. 
#all instances: if domain A appears twice in protein 4, count +2


import pandas as pd
import sys
import_path_base = '/storage1/hezscha/src/'
#import_path_base = '/home/henrike/Documents/PD_AS/src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
from PrismData import PrismParser, VariantData#VariantParser, VariantData
from Bio import AlignIO
import numpy as np
import glob
import os
import math

def read_from_prism(primsfile):
	parser = PrismParser()
	dataframe = parser.read(primsfile).dataframe
	meta_data = parser.read_header(primsfile)
	return meta_data, dataframe

def write_prism(metadata, dataframe, prism_file, comment=''):
	variant_dataset = VariantData(metadata, dataframe)
	parser = PrismParser()
	parser.write(prism_file, variant_dataset, comment_lines=comment)

prism_dir = "/storage1/shared/data/prism/"
outdir = '/storage1/hezscha/genome_proteome_map/results/domain_analysis/domain_lists/'

for d1 in os.scandir(prism_dir): 
	try:
		for d2 in os.scandir(d1.path):
			try:
				for d3 in os.scandir(d2.path):
					for prism_file in os.scandir(d3.path):
						if prism_file.name.endswith('_filled.txt') or prism_file.name.endswith('_SNPinv.txt'):
							continue
						
						elif prism_file.name.startswith('prism_uniprot_002_'):
							
							uniprot = prism_file.name.split('_')[-1].split('.')[0]
							print(uniprot)	
							metadata,df = read_from_prism(prism_file.path)
							
							if 'pfam_name' in df.columns:
								domains = df.pfam_name.unique()
								
								for dom in domains:
									#print(dom)
									#skip nan
									#figured out nan is not a string. I can't test with isnan because that fails when dom is a string 
									if not isinstance(dom, str):
										#print('Skipping', dom)
										continue
									
									list_file = os.path.join(outdir, dom+'.uniprot')
									#does a list already exist? if not, create it
									if not os.path.exists(list_file):
										with open(list_file, 'w') as OUT:
											print(uniprot, file = OUT)
									
									#else open in append mode
									else:
										with open(list_file, 'a') as OUT:
											print(uniprot, file = OUT)
		
			except NotADirectoryError:
				continue				
	except NotADirectoryError:
		continue						
