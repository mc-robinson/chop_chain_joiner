#!/usr/bin/env python
# File name: make_pdb_seq.py
# Author: Matt Robinson
# Date created: 5/23/2017
# Date last modified: 5/23/2017
# Python Version: 3.6
"""
Description:
This script extracts the sequence from a PDB file with missing residues. 
The output is a sequence file of the form [pdb_code].seq

Usage: make_pdb_seq.py pdb_filename

Note: the PDB file must be in the same directory as this script.
"""
import sys
import os
from modeller import *

#Get the PDB id from the file
pdb_id = os.path.splitext(os.path.basename(sys.argv[1]))[0]

# Get the sequence of the PDB file, and write to an alignment file
code = pdb_id

e = environ()
m = model(e, file=code)
aln = alignment(e)
aln.append_model(m, align_codes=code)
aln.write(file=code+'.seq')





