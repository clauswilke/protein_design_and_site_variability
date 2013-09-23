#!usr/local/bin/python
import sys, os, math, string, re, gzip, urllib, shutil, Bio, subprocess
import cStringIO 
from Bio.PDB import *
from Bio.PDB.DSSP import * 
from scipy.stats import pearsonr as pearson
import random as rnd
import numpy as np
from Bio.PDB.Polypeptide import *
from pylab import *
import matplotlib.pyplot as plt

#LAST UPDATED: June 24, 2013
#Description: This is a series of functions that are used in the Protein Design Project Analysis

EulerGamma = 0.57721566490153286060 #Euler Gamma Constant

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

#This is a dictionary that has the amino acid for the key and the max solvent accessibility for this amino acid
#THIS HAS BEEN UPDATED. I AM USING THE NEW THEORETICAL NUMBERS FROM THE 2013 AUSTIN, STEPHANIE, MATT, WILKE PAPER. 
residue_max_acc = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
                   'C': 158.0, 'Q': 223.0, 'E': 224.0, 'G': 104.0,  \
                   'H': 209.0, 'I': 197.0, 'L': 201.0, 'K': 237.0, \
                   'M': 218.0, 'F': 239.0, 'P': 159.0, 'S': 151.0, \
                   'T': 172.0, 'W': 282.0, 'Y': 263.0, 'V': 174.0}

identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076' :'2X1B','PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37','PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY','PF00355' :'1FQT','PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS','PF00542' :'1CTF','PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8','PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI','PF02823' :'1AQT','PF04002': '2QLC','PF07686' : '2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

#These three functions take an array and then create a formatted string that is used for printing data to a file.
def dump_csv_line(line):
    new_line = ""
    size = len(line)
    for x in range(0,size):
        new_line = new_line + str(line[x])
        if(x != size-1):
            new_line = new_line + ","

    new_line += "\n"
    return new_line
 
def dump_csv_line2(line):
    new_line = ""
    size = len(line)
    for x in range(0,size):
        new_line = new_line + str(line[x])
        if(x != size-1):
            new_line = new_line + " "

    new_line += "\n"
    return new_line

def dump_csv_line3(line):
    new_line = ""
    size = len(line)
    for x in range(0,size):
        new_line = new_line + str(line[x])
        if(x != size-1):
            new_line = new_line + "\t"

    new_line += "\n"
    return new_line

#This functions takes a file with a bunch of natural sequences that have been aligned and then returns a list of the sequences. 
def get_natural_sequences(file):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []

    fileparts = re.split("_",file)
    pdb_id = fileparts[0]
    chain_id = fileparts[1]
    
    file_data = open(file, "r") #This is were the sequences are stored.
    seq_data = file_data.readlines()
    file_data.close()
    string  = ''
    finished_sequence = ''
    for sequence in seq_data:
        if sequence[0] == '>':
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = ''
        else:       
            string = string + sequence.rstrip("\n")
    all_sequences.append(string)
    num_sequences = len(all_sequences)
    return all_sequences   

def get_cut_designed_sequences(file):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    fileparts = re.split("_",file)
    identity = fileparts[0]
    pdb_id = fileparts[1]
    chain_id = fileparts[2]
    index_filename = pdb_id + '_A.indices'
    file_data = open("/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/backrub_sequences/" + file, "r") #This is were the sequences are stored.
    seq_data = file_data.readlines()
    file_data.close()
    index_file = open("/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/structures/" + index_filename)
    index_data = index_file.readlines()
    index_file.close()
    indices = []
    mod_sequences = []
    for line in index_data:
        parts = re.split("\t", line)
        index = int(parts[2])
        indices.append(index)

    string  = ''
    finished_sequence = ''
    for sequence in seq_data:
        if sequence[0] == '>':
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = ''
        else:       
            string = string + sequence.rstrip("\n")
    all_sequences.append(string)
    num_sequences = len(all_sequences)

    for seq in all_sequences:
        new_seq = ''
        #print len(seq)
        #new_list = seq.split('')
        #num_acids = len(new_list)
        for i in indices:
            new_seq = new_seq + seq[i-1]
        mod_sequences.append(new_seq)
    return mod_sequences   

def get_cut_natural_sequences(sequence_identity, all_sequences):
    #all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    file = sequence_identity + '.align.80'
    #fileparts = re.split("_",file)
    #identity = fileparts[0]
    pdb_id = identity_dict[sequence_identity]
    chain_id = 'A'
    #pdb_id = fileparts[1]
    #chain_id = fileparts[2]
    index_filename = pdb_id + '_A.indices'
    file_data = open("/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/natural_alignments/" + file, "r") #This is were the sequences are stored.
    seq_data = file_data.readlines()
    file_data.close()
    index_file = open("/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/structures/" + index_filename)
    index_data = index_file.readlines()
    index_file.close()
    indices = []
    mod_sequences = []
    for line in index_data:
        parts = re.split("\t", line)
        index = int(parts[0])
        indices.append(index)
    '''
    string  = ''
    finished_sequence = ''
    for sequence in seq_data:
        if sequence[0] == '>':
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = ''
        else:       
            string = string + sequence.rstrip("\n")
    all_sequences.append(string)
    num_sequences = len(all_sequences)
    '''
    for seq in all_sequences:
        new_seq = ''
        #new_list = seq.split('')
        #num_acids = len(new_list)
        for i in indices:
            new_seq = new_seq + seq[i]
        mod_sequences.append(new_seq)
    return mod_sequences   

def get_cut_RSA_values(pdb_id, RSA_values, RSA_dict):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    #fileparts = re.split("_",file)
    #identity = fileparts[0]
    #pdb_id = fileparts[1]
    #chain_id = fileparts[2]
    index_filename = pdb_id + '_A.indices'
    #file_data = open("/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/backrub_sequences/" + file, "r") #This is were the sequences are stored.
    #seq_data = file_data.readlines()
    #file_data.close()
    #seq_data = file_data.readlines()
    #file_data.close()
    index_file = open("/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/structures/" + index_filename)
    index_data = index_file.readlines()
    index_file.close()
    indices = []
    mod_sequences = []
    for line in index_data:
        parts = re.split("\t", line)
        #index = int(parts[3])
        index = re.split(" ",parts[3])
        index = int(index[0])
        indices.append(index)

    #mod_RSA_values = []
    #RSA_indices = RSA_dict.keys()
    #for key in RSA_indices:
    #    mod_RSA_values.append(RSA_dict[key])
    #print RSA_indices
    mod_RSA_values = []
    #print len(RSA_values)
    #print indices
    for i in indices:
        #print i
        #print len(RSA_values)
        mod_RSA_values.append(RSA_dict[i])
    print "Length of cut RSA values: " + str(len(mod_RSA_values))
    return mod_RSA_values  



#This takes a file with a bunch of aligned natural and designed sequences and returns two lists with them seperated. 
def split_merged_sequences(file):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    fileparts = re.split("_",file)
    pdb_id = fileparts[0]
    chain_id = fileparts[1]
    temp = fileparts[2]
    file_data = open(file, "r")
    seq_data = file_data.readlines()
    file_data.close()
    string  = ''
    finished_sequence = ''
    for sequence in seq_data:
        #print sequence
        if sequence[0] == '>':
            #print sequence
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = ''
        else:       
            string = string + sequence.rstrip("\n")
    all_sequences.append(string)
    num_sequences = len(all_sequences)
    num_designed_seq = 500
    num_natural_seq = num_sequences - num_designed_seq
    for seq_index in xrange(0,num_natural_seq):
        natural_sequences.append(all_sequences[seq_index])
    for seq_index in xrange(num_natural_seq,num_sequences):
        designed_sequences.append(all_sequences[seq_index])
    return natural_sequences, designed_sequences, temp

#This takes a list of natural sequences and splits them in half. It returns two samples with the split sequences.
def split_natural_sequences(natural_sequences):
    seq_sample_1 = natural_sequences
    L = len(natural_sequences)
    L_sample = int(L/2)
    seq_sample_2 = rnd.sample(natural_sequences, L_sample)
    for seq in seq_sample_2: 
        seq_sample_1.remove(seq)
    return seq_sample_1,seq_sample_2

#Gets the number of times that a count value is seen in the array and returns a list with the count, countOfCounts data 
def get_AA_counts(distribution_data):
    count_data = []
    for x in distribution_data:
        if x != 0:
            num_appearances = distribution_data.count(x)
            count_data.append((x,num_appearances))
    count_data = list(set(count_data))
    return count_data

#This function takes a list of sequences and returns a list of lists that contain all of the frequency data for every site
def get_transformed_data(list_of_sequences):
    transformed_distribution = []
    RSA = []
    for seq in list_of_sequences:
        sequences = seq.strip()
    list_of_sequences.pop(0)

    new_data = []
    for line in list_of_sequences:
        element = line.split()
        element.pop(0)
        element.pop(0)
        element.pop(0)
        new_data.append(element)
    for data in new_data:
        new_elements = []
        for count in data:       
            new_count = int(count)
            new_elements.append(new_count)
        transformed_distribution.append(new_elements)
    return transformed_distribution

#This function takes the amino acid count data for the designed and corresponding natural sequences and returns a list of lists with all the AA data at each site. 
def get_AA_distribution(proteins):
    #print proteins
    protein_distribution = []
    transformed_protein_distribution = []
    num_AA = 0
    input = open(proteins, "r")
    protein_data = input.readlines()
    input.close()
    protein_distribution = get_transformed_data(protein_data)
    for site in protein_distribution:
        new_site = site
        num_AA = sum(new_site)
        aa_probs = []
        for count in new_site:   
            if count == 0:
                prob = float(1)/float(num_AA + 20)
            else:
                prob = float(count + 1)/float(num_AA + 20)
            aa_probs.append(prob)
        transformed_protein_distribution.append(aa_probs)
    return transformed_protein_distribution
 
#Returns a list of KL_Divergence values for two distributions. The inputs are two list of lists representing the two sequence alignments
def get_Kullback_Leibler(real_proteins, created_proteins):
    KL_Values = []
    KL_Number = 0
    real_array = array(real_proteins)
    created_array = array(created_proteins)
    created_num_residues, created_num_AA = created_array.shape
    num_residues,num_AA = real_array.shape
    for i in xrange(0, num_residues):
        real_values = real_proteins[i]
        created_values = created_proteins[i]
        KL_Number = 0
        for j in xrange(0,20):
            value = log(float(real_values[j])/float(created_values[j]))
            value = value*float(real_values[j])
            KL_Number = KL_Number + value
        KL_Values.append(KL_Number)
    return KL_Values
 
#This function calculates the mean of list of data points and returns the mean
def calculate_mean(data_points):
    mean_of_data = []
    new_data_points = []
    for point in data_points:
        new_point = float(point)
        new_data_points.append(new_point)
    
    sum_of_data = sum(new_data_points)
    num_elements = len(data_points)
    mean_of_data = float(sum_of_data)/float(num_elements)
    return mean_of_data
 
#This function takes a file and extracts the RSA values
def get_RSA_Values(protein_file):
    input = open(protein_file, "r")
    protein_data = input.readlines()
    input.close()
    RSA = []
    for seq in protein_data:
        sequence = seq.strip()
    protein_data.pop(0)
    new_data = []
    for line in protein_data:
        element = line.split()
        test = element.pop(0)
        test = element.pop(0)
        RSA_value = element.pop(0)
        RSA.append(RSA_value)
        new_data.append(element)
    return RSA

#The inputs to this functions are three lists with the RSA, entropy, and KL site data for a protein. This sorts the entropy, and 
#KL data by RSA into three categories: buried, partially buried, and surface. It then returns the data. 
def get_position_dependent_data(RSA_data, entropy_data, KL_data):
    num_sites = len(RSA_data)
    i = 0
    buried_KL_values = []
    buried_entropy_values = []
    buried_RSA_values = []

    intermediate_KL_values = []
    intermediate_entropy_values  = []
    intermediate_RSA_values = []

    surface_KL_values = []
    surface_entropy_values = []
    surface_RSA_values = []
    
    while i < num_sites:
        if (float(RSA_data[i])<0.05):
            buried_KL_values.append(float(KL_data[i]))
            buried_entropy_values.append(float(entropy_data[i]))
            buried_RSA_values.append(float(RSA_data[i]))
            i = i + 1
        elif (0.05<=float(RSA_data[i])<=0.25): 
            intermediate_KL_values.append(float(KL_data[i]))
            intermediate_entropy_values.append(float(entropy_data[i]))
            intermediate_RSA_values.append(float(RSA_data[i]))
            i = i + 1
        elif (float(RSA_data[i]) > 0.25):
            surface_KL_values.append(float(KL_data[i]))
            surface_entropy_values.append(float(entropy_data[i]))
            surface_RSA_values.append(float(RSA_data[i]))
            i = i + 1
        else:
            print "Problem in get_core_data"
            print "RSA value is: " + str(RSA_data[i])
            i = i + 1
    return buried_entropy_values, buried_KL_values, intermediate_entropy_values, intermediate_KL_values, surface_entropy_values, surface_KL_values 

#This takes a file with site amino acid count data and then calculates the entropy for sites. 
def get_native_entropy(protein_file):
    probs = get_AA_distribution(protein_file)
    entropy_values = []
    entropy_number = 0
    probs_array = array(probs)
    num_residues,num_AA = probs_array.shape
    
    for i in xrange(0, num_residues):
        probs_values = probs_array[i]
        prob_sum = sum(probs_values)
        entropy_number = 0
        for j in xrange(0,20):
            value = (float(probs_values[j])*log(float(probs_values[j])))
            entropy_number = entropy_number + value
        entropy_values.append(-entropy_number)
    return entropy_values

#This takes the file with the count data and calculates the entropy using a different entropy calculator. 
def get_entropy(protein_file):
    protein_distribution = []
    transformed_protein_distribution = []
    num_AA = 0
    input = open(protein_file, "r")
    protein_data = input.readlines()
    input.close()
    probs = get_transformed_data(protein_data)
    entropy_values = []
    entropy_number = 0
    probs_array = array(probs)
    new_probs = []
    new_probs_array = array(probs)
    for site in probs:
        new_site_array = []
        for count in site:
            new_count = count + 1
            new_site_array.append(new_count)
        new_probs.append(new_site_array)
    num_residues,num_AA = new_probs_array.shape
    count = 0
    for site in new_probs:
        frequencies = site
        entropy_number = Entropy_H_G(frequencies)
        entropy_values.append(entropy_number)
        if (entropy_number >2.998):
            print count
            print frequencies
            print "Entropy Number: " + str(entropy_number)  
        count = count + 1
    return entropy_values

#This takes a list of values that are string values and creates an array
def make_array(list_of_values):
    new_value_list = []
    for value in list_of_values:
        new_value = float(value)
        new_value_list.append(new_value)
    value_array = array(new_value_list)
    return value_array

#Opens a file with the slope and intercept data for all the proteins for the designed proteins. 
#It returns two dictionaries that map the slope and intercept data to the PDB names. 
def get_slope_intercept(file):
    slope_list = []
    intercept_list = []
    pdb_list = []
    input = open(file, "r")
    file_data = input.readlines()
    for line in file_data:
        data = re.split(",", line)
        pdb_line = data[0]
        PDB = re.split("_", pdb_line)
        PDB = PDB[4]
        pdb_list.append(PDB)
        intercept = data[1]
        slope = data[2]
        intercept_list.append(intercept)
        slope_list.append(slope)
    intercept_dict = dict(zip(pdb_list, intercept_list))
    slope_dict = dict(zip(pdb_list, slope_list))
    return intercept_dict, slope_dict

#Opens a file with the slope and intercept data for all the proteins for the natural proteins. 
#It returns two dictionaries that map the slope and intercept data to the PDB names. 
def get_slope_intercept_natural(file):
    slope_list = []
    intercept_list = []
    pdb_list = []
    input = open(file, "r")
    file_data = input.readlines()
    for line in file_data:
        data = re.split(",", line)
        pdb_line = data[0]
        PDB = re.split("_", pdb_line)
        PDB = PDB[5]
        pdb_list.append(PDB)
        intercept = data[1]
        slope = data[2]
        intercept_list.append(intercept)
        slope_list.append(slope)
    intercept_dict = dict(zip(pdb_list, intercept_list))
    slope_dict = dict(zip(pdb_list, slope_list))
    return intercept_dict, slope_dict
'''
#This calculates the RSA values for a PDB using DSSP. It returns a list of the Amino Acids and a list of their RSA values. 
def get_values(pdb_id, chain_id):
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
    pdbLocation = "/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/structures/" + searchPDB #This is the location of the pdb file that dssp will be parsing. 
    outputFile = pdb_id + "_" + chain_id + ".txt" 
    processString = 'dssp' + ' -i ' + '"' + pdbLocation +'"'  + ' -o pdbOutput.txt '
    process = subprocess.Popen(processString, shell = True, stdout = subprocess.PIPE)
    process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values
    input = open("pdbOutput.txt" , 'r')    
    fileContents = input.readlines()	
    string = fileContents[25]
    SAValue1 = string[13]
    SAList = [] #This is is the list which will store the SA values for each site
    AAList = [] #This is the list which will store the amino acid values for each site
    index = 0
    NoRSA = 0
    for line in fileContents:
        if index<25: #This skips the first few lines that do not have the SA value data
            index = index + 1
            continue
        else:  #Goes through each line with has the SA for each amino acid in order
            string = line #This stores the current line in the string "string"
            SAValue = string[35:39] #This stores the SA value for the current 
            AA = string[13] #This stores what the amino acid type is at the current  position
            number = int(SAValue)  #This turns the string fir the SA into an int type
            if AA !=( '!' or '*'): #This takes out the missing gaps that dssp might put in
                max_acc = residue_max_acc[AA] #This uses the dictionary to find the max SA for the amino acid at the current position (site)
                SAList.append(number/max_acc) #This divides the SA value for that position by the amx SA value position for tha amino acid. This normalizes the values and gives us the Relative Solvent Accessability (RSA) value. We this appends this value to the list of RSA values
                AAList.append(AA) #This appends the amino acid to the list
            else:
           	NoRSA = NoRSA + 1 #Counts the number of residues that did not have SA values
	    index = index + 1
    input.close() #Close the file with the dssp output file
    #os.remove('pdbOutput.txt') #Deletes the dssp output file 
    return (AAList, SAList) #Return the RSA values and the SAList
'''

#This calculates the RSA values for a PDB using DSSP. It returns a list of the Amino Acids and a list of their RSA values. 
def get_values(pdb_id, chain_id):
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
    pdbLocation = "/home/elj299/RSA_from_ceres/other/structs_all/" + searchPDB #This is the location of the pdb file that dssp will be parsing.
    structure = PDBParser().get_structure(pdb_id, pdbLocation)#Creates a structure object 
    model = structure[0]
    temp_file = 'new_pdb_tempfile.txt'
    Bio.PDB.Dice.extract(structure, chain_id, 0, 10000, temp_file)
    outputFile = pdb_id + "_" + chain_id + ".txt" 
    #processString = 'dssp' + ' -i ' + '"' + pdbLocation +'"'  + ' -o pdbOutput.txt '
    processString = 'dssp' + ' -i ' + '"' + temp_file +'"'  + ' -o pdbOutput.txt '
    process = subprocess.Popen(processString, shell = True, stdout = subprocess.PIPE)
    process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values
    input = open("pdbOutput.txt" , 'r')    
    fileContents = input.readlines()	
    string = fileContents[25]
    SAValue1 = string[13]
    SAList = [] #This is is the list which will store the SA values for each site
    AAList = [] #This is the list which will store the amino acid values for each site
    index = 0
    NoRSA = 0
    for line in fileContents:
        if index<25: #This skips the first few lines that do not have the SA value data
            index = index + 1
            continue
        else:  #Goes through each line with has the SA for each amino acid in order
            string = line #This stores the current line in the string "string"
            SAValue = string[35:39] #This stores the SA value for the current 
            AA = string[13] #This stores what the amino acid type is at the current  position
            number = int(SAValue)  #This turns the string fir the SA into an int type
            if AA !=( '!' or '*'): #This takes out the missing gaps that dssp might put in
                max_acc = residue_max_acc[AA] #This uses the dictionary to find the max SA for the amino acid at the current position (site)
                SAList.append(number/max_acc) #This divides the SA value for that position by the amx SA value position for tha amino acid. This normalizes the values and gives us the Relative Solvent Accessability (RSA) value. We this appends this value to the list of RSA values
                AAList.append(AA) #This appends the amino acid to the list
            else:
           	NoRSA = NoRSA + 1 #Counts the number of residues that did not have SA values
	    index = index + 1
    input.close() #Close the file with the dssp output file
    #os.remove('pdbOutput.txt') #Deletes the dssp output file 
    return (AAList, SAList) #Return the RSA values and the SAList


#This calculates the RSA values for a PDB using DSSP. It returns a list of the Amino Acids and a list of their RSA values. 
def get_noah_RSA_values(pdb_id, chain_id):
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
    pdbLocation = "/home/elj299/RSA_from_ceres/project/noah_sequence_analysis/structures/" + searchPDB #This is the location of the pdb file that dssp will be parsing. 
    structure = PDBParser().get_structure(pdb_id, pdbLocation)#Creates a structure object 
    model = structure[0]
    temp_file = 'new_pdb_tempfile.txt'
    Bio.PDB.Dice.extract(structure, chain_id, 0, 10000, temp_file)
    outputFile = pdb_id + "_" + chain_id + ".txt" 
    #processString = 'dssp' + ' -i ' + '"' + pdbLocation +'"'  + ' -o pdbOutput.txt '
    processString = 'dssp' + ' -i ' + '"' + temp_file +'"'  + ' -o pdbOutput.txt '
    process = subprocess.Popen(processString, shell = True, stdout = subprocess.PIPE)
    process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values
    input = open("pdbOutput.txt" , 'r')    
    fileContents = input.readlines()	
    string = fileContents[25]
    SAValue1 = string[13]
    SAList = [] #This is is the list which will store the SA values for each site
    AAList = [] #This is the list which will store the amino acid values for each site
    residue_list = [] #This is a list that stores the residue postion numbers
    index = 0
    NoRSA = 0
    for line in fileContents:
        if index<25: #This skips the first few lines that do not have the SA value data
            index = index + 1
            continue
        else:  #Goes through each line with has the SA for each amino acid in order
            string = line #This stores the current line in the string "string"
            SAValue = string[35:39] #This stores the SA value for the current 
            res_pos = string[6:10]
            AA = string[13] #This stores what the amino acid type is at the current  position
            number = int(SAValue)  #This turns the string for the SA into an int type
            if AA !=( '!' or '*'): #This takes out the missing gaps that dssp might put in
                max_acc = residue_max_acc[AA] #This uses the dictionary to find the max SA for the amino acid at the current position (site)
                SAList.append(number/max_acc) #This divides the SA value for that position by the amx SA value position for tha amino acid. This normalizes the values and gives us the Relative Solvent Accessability (RSA) value. We this appends this value to the list of RSA values
                residue = int(res_pos)
                AAList.append(AA) #This appends the amino acid to the list
                residue_list.append(residue)
            else:
           	NoRSA = NoRSA + 1 #Counts the number of residues that did not have SA values
	    index = index + 1
    RSA_dict = dict(zip(residue_list,SAList))
    input.close() #Close the file with the dssp output file
    #os.remove('pdbOutput.txt') #Deletes the dssp output file 
    return (AAList, SAList, RSA_dict) #Return the RSA values and the SAList

#This takes a file with a bunch of generated quantities (mean RSA, mean KL value, mean entropy, correlation between RSA and Entropy) 
#for a protein and returns lists with this data
def get_mean_designed_data(file_of_data):
    designed_mean_RSA_values = []
    designed_mean_entropy_values = []
    designed_cor_entropy_RSA_values = []
    designed_mean_KL_values = []
    designed_data = []
    protein_file = open(file_of_data, "r")
    designed_protein_data = protein_file.readlines()
    protein_file.close()
    header = designed_protein_data.pop(0)

    for line in designed_protein_data:  
        data = re.split("\t", line)
        designed_data.append(data)
    
    for data in designed_data:
        designed_mean_RSA_values.append(data[2])
        designed_mean_entropy_values.append(data[3])
        designed_cor_entropy_RSA_values.append(data[4])
        designed_mean_KL_values.append(data[5])
    return designed_mean_RSA_values, designed_mean_entropy_values, designed_cor_entropy_RSA_values, designed_mean_KL_values

#This takes a file with a bunch of generated quantities (RSA, KL value, Entropy) for a protein and returns lists with this data
def get_designed_graph_data(generated_data_file):
  designed_file = open(generated_data_file)
  data = designed_file.readlines()
  designed_data = []
  designed_file.close()
  for line in data:
      file_data = re.split("\t",line)
      designed_data.append(file_data)

  designed_entropy = designed_data[0]
  designed_entropy.pop(0)
  designed_entropy_array = make_array(designed_entropy)

  designed_RSA = designed_data[1]
  designed_RSA.pop(0)
  designed_RSA_array = make_array(designed_RSA)

  designed_KL = designed_data[2]
  designed_KL.pop(0)
  designed_KL_array = make_array(designed_KL)   
  return designed_entropy,designed_RSA,designed_KL    

def get_mean_ordered_designed_data(file_of_data):
    designed_mean_RSA_values = []
    designed_mean_entropy_values = []
    designed_cor_entropy_RSA_values = []
    designed_mean_KL_values = []
    designed_intercept_values = []
    designed_slope_values = []
    designed_data = []
    protein_file = open(file_of_data, "r")
    designed_protein_data = protein_file.readlines()
    protein_file.close()
    header = designed_protein_data.pop(0)

    for line in designed_protein_data:  
        data = re.split("\t", line)
        designed_data.append(data)
    
    for data in designed_data:
        designed_mean_RSA_values.append(data[2])
        designed_mean_entropy_values.append(data[3])
        designed_cor_entropy_RSA_values.append(data[4])
        designed_mean_KL_values.append(data[5])
        designed_intercept_values.append(data[6])
        designed_slope_values.append(data[7])

    return designed_mean_RSA_values, designed_mean_entropy_values, designed_cor_entropy_RSA_values, designed_mean_KL_values, designed_intercept_values, designed_slope_values

#This is a G() calculator that is used to calculate entropy
def Gi(n): #Function Written By Claus Wilke
    '''Helper function needed for entropy estimation.
    Defined by Grassberger 2003. http:/arxiv.org/abs/physics/0307138
    '''  
    if n == 0:
        return 0
    if n == 1:
        return -EulerGamma - math.log(2)
    if n == 2:
        return 2-EulerGamma - math.log(2)
    if (n % 2) == 1:
        return Gi( n-1 )
    return Gi(n-2) + 2./(n-1)

#This is used in the entropy calculator
def Entropy_H_G(list_of_frequency_counts): #Function Written By Claus Wilke
    '''Best entropy estimator according to Grassberget 2003,
    http:/arxiv.org/abs/physics/0307138
    '''  
    z = list_of_frequency_counts
    N = sum(z) # total number of observations
    return math.log(N) - (1./N)* sum([n*Gi(n) for n in z])
