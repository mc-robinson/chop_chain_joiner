#!/usr/bin/env python
# File name: parse_chop.py
# Author: Matt Robinson
# Date created: 5/26/2017
# Date last modified: 5/26/2017
# Python Version: 3.6
"""
Description:

Usage: python chop_make_alignment.py chop.log pdbfile.pdb pdbseq.seq fastafile.fasta

Note: the PDB file must be in the same directory as this script.
"""
import sys
import re
import os 
from biopandas.pdb import PandasPdb
import pandas as pd

with open(sys.argv[1]) as chop_log:
	chop_data = chop_log.readlines()

with open(sys.argv[3]) as pdb_seq_file:
	pdb_seq_data = pdb_seq_file.readlines()

with open(sys.argv[4]) as fasta_file:
	full_seq_data = fasta_file.readlines()

ppdb = PandasPdb()
ppdb.read_pdb(sys.argv[2]) #need to fill this in.

pdb_id = os.path.splitext(os.path.basename(sys.argv[2]))[0]


def main():
	#get the lists of disconnections
	missing_atom_ls = get_atoms(chop_data)

	pdb_seq = get_pdb_seq(pdb_seq_data)
	full_seq = get_full_seq(full_seq_data)

	#now trim the full sequence so that the non-crystallized ends are removed
	#full_seq = trim_full_seq(pdb_seq, full_seq)
	print(pdb_seq)
	print('\n')
	print(full_seq)

	#need to trim off ligand of pdb_seq
	pdb_seq = trim_ligand_seq(pdb_seq)[0]
	print(pdb_seq)
	ligand_str = trim_ligand_seq(pdb_seq)[1]
	print(ligand_str)

	#break the two seqs up chain by chain
	pdb_seq_l = break_into_chains(pdb_seq)
	#print(pdb_seq_l)
	full_seq_l = break_into_chains(full_seq)
	#print(full_seq_l)

	#check if they have the same number of chains
	if not(len(pdb_seq_l) == len(full_seq_l)):
			print("Fasta sequence has incorrect number of chains")
			sys.exit(2)

	#loop over the disconnections
	for l in missing_atom_ls:
		seq_l = get_missing_seq(l)[0] #the seq residues flanking the break
		gap_res_l = get_missing_seq(l)[1] #residues b/w which break occurs
		seq_before = get_missing_seq(l)[2] #seq before break
		seq_after = get_missing_seq(l)[3] #seq after break

		#make into one letter seq to work with fasta file
		seq = make_one_letter(seq_l)
		gap_res = make_one_letter(gap_res_l)
		seq_before = make_one_letter(seq_before)
		seq_after = make_one_letter(seq_after)

		print(seq)
		print(gap_res)
		print(seq_before)
		print(seq_after)

		#loop over the chains
		for i in range(len(pdb_seq_l)):
			#use seqs defined above to test if in this chain
			test_seq = seq
			test_seq_before = seq_before
			test_seq_after = seq_after

			#trim the full sequence so that the non-crystallized ends are removed
			#also check if it is last sequence and need to take off '* character'
			if (not(i == len(pdb_seq_l) -1)):
				full_seq_l[i] = trim_full_seq(pdb_seq_l[i], full_seq_l[i])
			else:
				full_seq_l[i] = trim_full_seq(pdb_seq_l[i][0:-1], full_seq_l[i][0:-1])
				full_seq_l[i] = full_seq_l[i] + '*' #add back the '*' end character

			#if there are hetatoms present, need to deal with that
			if ('.' in pdb_seq_l[i]):
				test_seq, test_seq_before, test_seq_after = get_hetatm_seq(seq, pdb_seq_l[i], gap_res)

			# check that this is the chain with the disconnection
			if (test_seq in pdb_seq_l[i]):
				print('hi')
				#now add the gaps in the string
				pdb_seq_l[i] = add_gaps(test_seq_before,test_seq_after,full_seq_l[i],pdb_seq_l[i])

	# put seqs back together
	pdb_seq = '/'.join(pdb_seq_l)
	full_seq = '/'.join(full_seq_l)

	#add last character to full_seq
	full_seq = full_seq

	#split both of the seqs into 75 character substrings for .pir format
	pdb_seq_substr = split_sequence(pdb_seq)
	full_seq_substr = split_sequence(full_seq)

	#write the output file into the local directory
	aln_file = open("alignment.ali","w+")

	for line in pdb_seq_data[0:3]:
		aln_file.write(line)

	for substr in pdb_seq_substr:
		aln_file.write(substr + '\n')

	aln_file.write('>P1;' + pdb_id + '_fill' + '\n') #header for full seq
	aln_file.write('sequence:::::::::'+'\n')

	for substr in full_seq_substr[0:-1]:
		aln_file.write(substr + '\n')
	aln_file.write(full_seq_substr[-1]) #don't want to end with newline character

	aln_file.close()


#find the atoms between which residues are missing.
def get_atoms(chop_data):
	atom_ls = []
	for line in chop_data:
		if (line[0:16] == '   disconnection'):
			atoms = []
			atoms_str = line.split('atoms ',1)[1]
			atoms_str = atoms_str.split(' is',1)[0]
			atom1 = atoms_str.split(', ',1)[0]
			atom2 = atoms_str.split(', ',1)[1]
			atoms.append(int(atom1))
			atoms.append(int(atom2))
			atom_ls.append(atoms)
	return atom_ls

def get_pdb_seq(pdb_seq_data):
	pdb_seq = ""
	seq_data = pdb_seq_data[3:]
	for line in seq_data:
		line = line.rstrip()
		pdb_seq = pdb_seq + line
	#remove the break character '/'
	#pdb_seq = re.sub('/','',pdb_seq)
	return pdb_seq

def get_full_seq(fasta_data):
	full_seq = ''
	for line in fasta_data:
		if not(line[0] == '>'):
			line = line.rstrip() #I could have just messed up some crap with this
			full_seq = full_seq + line
		if (line[0] == '>'): #if have new chain
			full_seq = full_seq + '/' #to indicate chain break
	full_seq = full_seq[1:] #get rid of beginning '/'
	return full_seq + '*'

def trim_full_seq(pdb_seq, full_seq):
	pdb_letters = pdb_seq
	#trim off the front
	for i in range(len(full_seq)):
		if (full_seq[i:i+5] == pdb_letters[0:5]):
			full_seq = full_seq[i:]
			break #needed so proper length
	#trim off the back (make sure * is not present or will mess it up)
	for j in reversed(range(len(full_seq)+1)):
		if (full_seq[j-5:j] == pdb_letters[(len(pdb_letters)-5):]):
			full_seq = full_seq[0:j]
			break
	return full_seq

def get_missing_seq(missing_atoms_l):
	#make biopandas dataframe from pdb
	#ppdb.read_pdb(sys.argv[2]) #need to fill this in.

	l = missing_atoms_l # for convenience of writing

	# when selecting residues, need to know if at end of chain
	if (l[0]-40 < 0):
		lower_bound = 0
	else:
		lower_bound = l[0]-40

	if (l[1]+40 >= ppdb.df['ATOM'].shape[0]):
		upper_bound = ppdb.df['ATOM'].shape[0] - 1
	else:
		upper_bound = l[1]+40

	# get lists of the residues (and their numbers) before and after the missing residues
	res_before_l = list(ppdb.df['ATOM'].loc[lower_bound:l[0]]['residue_name'])
	resnum_before_l = list(ppdb.df['ATOM'].loc[lower_bound:l[0]]['residue_number'])

	res_after_l = list(ppdb.df['ATOM'].loc[l[1]:upper_bound]['residue_name'])
	resnum_after_l = list(ppdb.df['ATOM'].loc[l[1]:upper_bound]['residue_number'])

	# add the residue names to the list, and make sure only once per residue
	res_before = []
	res_before.append(res_before_l[0])
	for i in range(1,len(res_before_l)):
		if ((not(resnum_before_l[i] == resnum_before_l[i-1]))): ##is this wrong
			res_before.append(res_before_l[i])

	# add the residue names to the list, and make sure only once per residue
	res_after = []
	res_after.append(res_after_l[0])
	for i in range(1,len(res_after_l)):
		if ((not(resnum_after_l[i] == resnum_after_l[i-1]))):
			res_after.append(res_after_l[i])


	merged_seq = res_before + res_after
	gap_res = [res_before[-1],res_after[0]]
	return merged_seq, gap_res, res_before, res_after

def make_one_letter(res_l):
	# Amino acid 3-to-1 letter code dictionary.
	aaDict = dict({'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'})
	#make the lift one letter
	one_letter_code_seq = []
	for res in res_l:
		one_letter_code_seq.append(aaDict[res])

	#now make the list one string
	seq_str = ''.join(one_letter_code_seq)
	return seq_str

def get_hetatm_seq(merged_seq,pdb_seq,gap_res):
	for i in range(len(merged_seq), len(pdb_seq)):
		test_seq = pdb_seq[0:i].replace('.','')
		if (test_seq[len(test_seq)-len(merged_seq):] == merged_seq):
			merged_seq = pdb_seq[i-len(merged_seq):i]
			print(i)
			break

	seq_before = merged_seq.split(gap_res)[0] + gap_res[0]
	seq_after = gap_res[1] + merged_seq.split(gap_res)[1]

	return merged_seq, seq_before, seq_after



#def find_missing_length(seq_before,full_seq):
	#find where to insert gap
	#gap_idx = full_seq.find(seq_before)

def insert_dashes(string, index, num):
	dash_str = ''
	for i in range(num):
		dash_str = dash_str + '-'
	return string[:index] + dash_str + string[index:]

#Add dashes where the missing residues are in the pdb_seq
def add_gaps(seq_before,seq_after,full_seq,pdb_seq):
	#find where gap begins, need to add back length because gives idx at beginning of str
	begin_gap_idx = full_seq.find(seq_before) + len(seq_before) #probably wrong
	# find where gap ends
	end_gap_idx = full_seq.find(seq_after)
	#compute number of dashed
	num_dashes = end_gap_idx - begin_gap_idx
	#find where to insert gap in pdb sequence
	insert_idx = pdb_seq.find(seq_before + seq_after) + len(seq_before)
	#make new pdb_seq
	pdb_seq = insert_dashes(pdb_seq,insert_idx,num_dashes)

	return pdb_seq

def split_sequence(seq):
	return re.findall('.{1,75}', seq)

def break_into_chains(seq):
	return seq.split('/')

def trim_ligand_seq(hetatm_seq):
	ligand_str = ''
	seq = hetatm_seq[0:-1] #take off the final *
	new_seq = seq
	for j in reversed(range(len(seq))):
		if not(seq[j] == '/' or seq[j] == '.'):
			new_seq = seq[0:j+1]
			ligand_str = seq[j+1:]
			break
	new_seq = new_seq + '*'
	return new_seq, ligand_str


if __name__ == "__main__":
	main()
