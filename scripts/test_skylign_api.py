
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
import requests

#parse arguments
################################################################################
#https://stackoverflow.com/questions/50021282/python-argparse-how-can-i-add-text-to-the-default-help-message
parser = argparse.ArgumentParser()
parser.add_argument('-dom', dest="dom", required = True, help='The domain to run for.')
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")
args = parser.parse_args()
################################################################################


dom_name = args.dom
base_dir = '/storage1/hezscha/genome_proteome_map/domain_variant_analysis/results/'

aln_file = os.path.join(base_dir, dom_name, 'alignment', 'aln_seed_instance') 
print(aln_file)

server = 'http://skylign.org'
r = requests.get(server, headers={ "Content-Type" : "application/json"}, files={ 'file': aln_file,'processing': 'hmm_all'})
print(r)
decoded = r.json()
print(decoded)
print(decoded.keys)
