#!usr/local/bin/python
import sys, os, math, string, re, gzip, urllib, shutil
import cStringIO 
import subprocess
 
searchStr = "[a-zA-Z0-9_\.\-]*" + ".fasta"
#find all csv files that match the search string
for path, names, filename in os.walk('.',False):
	for file in filename:
		if(re.search(searchStr, file)!=None):
			fileparts = re.split("_",file)
			pdb_id = fileparts[0].upper()
			#print pdb_id
			chain_id = fileparts[1]
			chain_id = chain_id[0]
			print "Processsing file: " + file	
			print "PDB: " + pdb_id
			print "CHAIN: " + chain_id
				  
			#fp = open(file,"r")
			#data = fp.readlines()
			#fp.close()
			outputName = pdb_id + '_' + chain_id + '_Aligned_Sequences.fasta'
			#print "Something" + pdb_id
			fileLocation = "/Users/Eleisha/Desktop/FastaFiles/" + file 
			alignedSeqLocation = "/Users/Eleisha/Desktop/FastaFiles/alignedSequences/" + outputName
			#print pdbLocation
			processString = 'mafft ' + fileLocation + " > " + alignedSeqLocation
			print processString
			process = subprocess.Popen(processString, shell = True, stdout = subprocess.PIPE)
			process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values