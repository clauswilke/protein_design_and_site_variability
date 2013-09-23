#!usr/local/bin/python
import sys, os, math, string, re, gzip, urllib, shutil, Bio
import cStringIO 
from Bio.PDB import *
from Bio.PDB.DSSP import *
from numpy import *   

#UPDATED: January 28 2012
#Description: This is a code that takes the unaligned_fasta files from the duncan analysis, uses the PDB File to find the original sequence and then creates a new file PDB_chain_alignedSequences.fasta. Note: although called aligned, these sequences are not technically aligned. (I know. This needs to be changed.)


resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',         
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }        

pdb_id = ''
chain_id = ''
badProteins = []
#This string that searched to find all files in the directory.
searchStr = "[a-zA-Z0-9_\.\-]*" + "_unaligned_sequences" +  ".fasta"

print searchStr
print "/n"

debug = 0

#find all csv files that match the search string
for path, names, filename in os.walk('.',False):
    for file in filename:
        if(re.search(searchStr, file)!=None):
			
            fileparts = re.split("_",file)
            pdb_id = fileparts[0].upper() #Store PDB Name
            pdb_id = pdb_id[0:4]
            chain_id = fileparts[0].upper()
            chain_id = chain_id[4]  #Store chain
	
            print "Processsing file: " + file	
            print "PDB: " + pdb_id
            print "CHAIN: " + chain_id
			
            fp = open(file,"r")
            data = fp.readlines()
            fp.close()		
				
            fileName = pdb_id + "_unaligned_sequences" + ".fasta"
            
            #Creates a filename for the output file 
            if chain_id == '-': 
                sequenceHeader = ">" + pdb_id + "_A" + "\n"
                newSeqFilename =  "/home/elj299/duncan_sequences/alignedSequences/" + pdb_id + "_A" + ".fasta"
            else:
                sequenceHeader = ">" + pdb_id + "_" + chain_id + "\n"
                newSeqFilename = "/home/elj299/duncan_sequences/alignedSequences/" + pdb_id + "_" + chain_id + ".fasta"
                
                #This is the code block that gets the sequence
                searchPDB = pdb_id + ".pdb"  
                pdbLocation = "/home/elj299/duncan_sequences/duncan_structs/" + searchPDB #This is the location of where the PDBS are location on your machine
                structure = PDBParser().get_structure(pdb_id, pdbLocation)#Creates a structure object 
                model = structure[0]
                if chain_id != '-':
                    chain = model[chain_id]
                    print chain
                else:
                    chain = model['A']
                residue_list = []

                for r in chain:
                    residue = r.get_resname()
                    #print residue
                    residue_list.append(residue)
                    #Gets all of the residues in the protein
                print len(residue_list)

                sequence = []
                sequenceString = ''
                for AA in residue_list:
                    if(AA in resdict):
                        aminoAcid = resdict[AA] #This changes the three-letter code to a one-letter code sequence
                        sequence.append(aminoAcid)
                        sequenceString = sequenceString + aminoAcid
                    else:
                        if (pdb_id in badProteins):
                            continue
                        else:
                            badProteins.append(pdb_id)
                        #This code blocks executes if there is a weird amino acid in the file that is not one of regular ones.
                out = open(newSeqFilename, "w")
                out.write(sequenceHeader)
                out.write(sequenceString + "\n\n")
                for sequence in data:
                    out.write(sequence)
                out.close

print badProteins
#This is a list of all the proteins that for some reason had weird sequences. 
badPDBS = open("Bad_Proteins.txt", "w")
for protein in badProteins:
    badPDBS.write(protein + "\n")
badPDBS.close()
