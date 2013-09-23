#!usr/local/bin/python
import sys, os, math, string, re, gzip, urllib, shutil, Bio, subprocess
import cStringIO 
import numpy as np
import analysis_functions as af
import random 

#Date Last Modified:  June 19, 2013
#Description: This file takes the designed and natural sequences for each protein at each temperature, calculates the frequency data and prints that info to some files. 

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',                   
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',                   
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',                   
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }          

rows = ["\"site\"", "\"RSA\"", "\"aa1\"", "\"aa2\"", "\"aa3\"", "\"aa4\"", "\"aa5\"", "\"aa6\"", "\"aa7\"", "\"aa8\"", "\"aa9\"", "\"aa10\"", "\"aa11\"", "\"aa12\"", "\"aa13\"", "\"aa14\"", "\"aa15\"", "\"aa16\"", "\"aa17\"", "\"aa18\"", "\"aa19\"", "\"aa20\""]   
 
identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076' :'2X1B','PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37','PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY','PF00355' :'1FQT','PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS','PF00542' :'1CTF','PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8','PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI','PF02823' :'1AQT','PF04002': '2QLC','PF07686' : '2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

#These are the absolute PATH to the UCSF sequence dataset files (structures, designed sequences, and natural sequences)
noah_designed_sequence_path = "/Users/Eleisha/Documents/Wilke_Lab/Project_1/project/sequences/noah_sequences/backrub_sequences/"
noah_natural_sequence_path = "/Users/Eleisha/Documents/Wilke_Lab/Project_1/project/sequences/noah_sequences/natural_alignments/"
noah_structure_path = "/Users/Eleisha/Documents/Wilke_Lab/Project_1/project/structures/noah_structures/"

temps = ["0.0", "0.3","0.6", "0.9", "1.2", "1.8", "2.4"]
chain_id = 'A'
identities = identity_dict.keys()
mean_designed_protein_length = 0
designed_protein_lengths = []
for temp in temps:
    for identity in identities:
        pdb_id = identity_dict[identity]
        file = identity + '.align.80'
        designed_filename = identity + "_" + pdb_id + "_" + chain_id + "_" + temp + ".fasta"
        designed_file = designed_filename
        fileparts = re.split("_",file)
        print "Processsing file: " + file	
        print "PDB: " + pdb_id
        print "CHAIN: " + chain_id
        searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
        pdbLocation = noah_structure_path + searchPDB #This is the location of the pdb file
        natural_sequences =  af.get_natural_sequences_noah(file) #Gets a list with the natural sequences 
        cut_natural_sequences = af.get_cut_natural_sequences(identity, natural_sequences)
        print "Length of Cut Natural Sequences: " + str(len(cut_natural_sequences[0]))
        ancestor = cut_natural_sequences[0] #Gets the "ancestral sequence"
        ancestor_length = len(ancestor) #Grab the length of the ancestral sequence
        print "Ancestor Length: " + str(len(natural_sequences[0]))
        counter = 0
        gaps = 0
        gap_locations = []

        #Counts the gaps within the ancestral sequence in the alignment. We must take the gaps out before counting the amino acids
        while (counter < ancestor_length):
            acid = ancestor[counter]
            if acid == '-': 
                gaps = gaps + 1
                gap_locations.append(counter) #This is an array that tracks which residues have gaps
            counter = counter + 1

        #This section takes out the ancestral gaps from all the aligned sequences and then writes the files
        #to the results_PDB_ID_CHAIN_ID.csv 
        natural_pdb_file_title = "results_natural_" + pdb_id + "_" + chain_id + ".csv"
        natural_out_sequences = open(natural_pdb_file_title,"w") 

        #This takes all of the gaps out of the alignment
        for sequence in cut_natural_sequences: #For each sequence in the the file with the sequences 
            new_sequence = '' 
            index = 0
            while (index < ancestor_length): #Take out all of the gaps within the aligned proteins
                if index in gap_locations:
                    index = index + 1
                    continue
                else:
                    new_sequence = new_sequence + sequence[index]
                    index = index + 1
            natural_out_sequences.write(new_sequence + "\n")
        natural_out_sequences.close()

        #list of files that contain L vs RSA data
        data = []

        #list of sequences with missing RSA values
        bad_list = []

        #Get all the RSA values using DSSP
        seq_data = af.get_noah_RSA_values(pdb_id, chain_id)
        #print seq_data
        RSAValues = []	
        RSA = seq_data[1]
        RSA_dict = seq_data[2]
        #print "RSA Length from dssp script: " + str(len(RSA))
        new_RSA = af.get_cut_RSA_values(pdb_id, RSA, RSA_dict)
        
        index = 0

        natural_data_arr = cut_natural_sequences       
        designed_data_arr = af.get_cut_designed_sequences(designed_file)
        #RSA = cut_RSA(natural
        #filepointer for the partial result
        fpW_natural = open("results_array_natural_" + pdb_id + "_" + chain_id + ".csv","w")
        fpW_designed = open("results_array_" + pdb_id + "_" + chain_id + "_" + temp + ".csv","w")

        #write the RSA values:
        fpW_natural.write(af.dump_csv_line(new_RSA))
        fpW_designed.write(af.dump_csv_line(new_RSA))

        if(len(new_RSA) != len(natural_data_arr[0].strip())):
            print "Error !!!!!!!!!"
            bad_list.append(natural_pdb_file_title)

        if(len(new_RSA) != len(natural_data_arr[0].strip())):
            print "Error !!!!!!!!!"
            print "Length of rsa values {:0}".format(len(RSA))
            print "Length of sequence {:0}".format(len(natural_data_arr[0].strip()))
            print "Length of designed sequence {:0}".format(len(designed_data_arr[0].strip()))

        #Re-format sequences
        for line in natural_data_arr:
            fpW_natural.write(af.dump_csv_line(line.strip()))      
        fpW_natural.close()

        for line in designed_data_arr:
            fpW_designed.write(af.dump_csv_line(line.strip()))      
        fpW_designed.close()

        #search string to use
        file = "results_array_natural_" + pdb_id + "_" + chain_id + ".csv" 
        print file 

        debug = 0

        #grab the file
        fp = open(file,"r")
        natural_data = fp.readlines()
        fp.close()

        #search string to use
        file = "results_array_" + pdb_id + "_" + chain_id +"_" + temp + ".csv" 
        print file 

        debug = 0

        #find all csv files that match the search string

        #grab the file
        fp = open(file,"r")
        designed_data = fp.readlines()
        fp.close()

        #split the dataset by RSA
        natural_size = len(natural_data)
        for i in range(0,natural_size):
            natural_data[i] = re.split(",",natural_data[i].strip())
            #print natural_data[i]
        if(debug):
            print len(natural_data[i])
        natural_data = np.array(natural_data, dtype = object)
        #print natural_data
        #Make arrays full of all the sequences from the files
        designed_size = len(designed_data)
        for i in range(0,designed_size):
            designed_data[i] = re.split(",",designed_data[i].strip())

        if(debug):
            print len(designed_data[i])
        designed_data = np.array(designed_data, dtype = object)

        #snag the RSA values
        RSA = natural_data[0]

        #get a list of char amino acid codes
        AA = resdict.values()
        #print "Dimensions of array {}".format(natural_data.ndim)
        
        n,m = natural_data.shape
        n_d, m_d = designed_data.shape 
        #print "Natural Protein Length", m
        #print "Designed Protein Length", m_d
        designed_protein_lengths.append(m_d)
        
        #Open .dat files that you will write the results to.
        natural_fpW = open("align_natural_data_array_" + pdb_id + "_" + chain_id + ".dat","w")
        natural_fpW.write(af.dump_csv_line2(rows))
        designed_fpW = open("align_data_array_" + pdb_id + "_" + chain_id + "_" + temp + ".dat","w")
        designed_fpW.write(af.dump_csv_line2(rows))

        counter = 0
        natural_aaSum = 0
        designed_aaSum = 0
        aaCutOff = 0 #This is cut-off. If a site has less than this percent of amino acids left after taking out the gaps do not use that site in the analysis. 			  
        for i in range(0,m): #This blocks calculates the site frequency data 
            natural_aaCount = [] 
            designed_aaCount = []
            #covert back to list so we can use the "count" method
            try:
                natural_aaList = list(natural_data[1:,i])
                designed_aaList = list(designed_data[1:,i])
            except IndexError:
                print "--------------"
                print "Fudge"
                print i
                print data.ndim
                print data.shape
                print "--------------"
                quit()

            for aa in AA:
                try:
                    natural_aaCount.append(natural_aaList.count(aa))
                    designed_aaCount.append(designed_aaList.count(aa))
                except ValueError,IndexError:
                    print "-----------------------"
                    print "Blargtastic!"
                    print aaList
                    print "-----------------------"
                    quit()

            #Print the formatted frequency count data to a file
            natural_outStr = "\"" + str(counter) + "\" "+ str(counter) + " " + str(RSA[i]) + " " + af.dump_csv_line2(natural_aaCount)  
            designed_outStr = "\"" + str(counter) + "\" "+ str(counter) + " " + str(RSA[i]) + " " + af.dump_csv_line2(designed_aaCount)
            counter = counter + 1
            natural_fpW.write(natural_outStr)
            designed_fpW.write(designed_outStr)

        print "Processed: %s\n" % file
        natural_fpW.close()
        designed_fpW.close() 

print designed_protein_lengths
mean_designed_protein_length = np.mean(designed_protein_lengths)
print "Noah Mean_Designed_Length: ", mean_designed_protein_length