import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from scipy.stats import pearsonr
from scipy.stats import ttest_rel
import analysis_functions


#Last Date Updated: August 7, 2013
#Description: This is a script that takes all of the *.dat files with the generated AA count and RSA values and then calculates the necessary values for comparision. Ex. Entropy, KL Divergence

#Get files 
#This searches for all of the unaligned sequences 
searchStr = "align_natural_data_array_" + "[a-zA-Z0-9_\.\-]*" +  ".dat"
#find all csv files that match the search string
temp = 0.0
temp_array = [0.0, 0.1, 0.3, 0.6, 0.9, 1.2]
data = []
count = 1 
#These array that are used to store various types of data
natural_mean_RSA_values = []
natural_mean_entropy_values = []
natural_cor_entropy_RSA_values = []

#Mean RSA
designed_mean_RSA_values_00 = []
designed_mean_RSA_values_01 = []
designed_mean_RSA_values_03 = []
designed_mean_RSA_values_06 = []
designed_mean_RSA_values_09 = []
designed_mean_RSA_values_12 = []
designed_mean_RSA_values_003 = []

natural_mean_RSA_values = []
natural_mean_entropy_values = []

#KL-Divergence
KL_array_00 = []
KL_array_01 = []
KL_array_03 = []
KL_array_06 = []
KL_array_09 = []
KL_array_12 = []
KL_array_003 = []

#Mean KL-Divergence
designed_mean_KL_values_00 = []
designed_mean_KL_values_01 = []
designed_mean_KL_values_03 = []
designed_mean_KL_values_06 = []
designed_mean_KL_values_09 = []
designed_mean_KL_values_12 = []
designed_mean_KL_values_003 = []

#Mean Entropy
designed_mean_entropy_values_00 = []
designed_mean_entropy_values_01 = []
designed_mean_entropy_values_03 = []
designed_mean_entropy_values_06 = []
designed_mean_entropy_values_09 = []
designed_mean_entropy_values_12 = []
designed_mean_entropy_values_003 = []

#The correlation between RSA and Entropy at sites
designed_cor_entropy_RSA_values_00 = [] 
designed_cor_entropy_RSA_values_01 = []
designed_cor_entropy_RSA_values_03 = []
designed_cor_entropy_RSA_values_06 = []
designed_cor_entropy_RSA_values_09 = []
designed_cor_entropy_RSA_values_12 = []
designed_cor_entropy_RSA_values_003 = []

split_KL_array = []
natural_mean_split_KL_values = []

#Mean KL- Divergence for Buried Sites
buried_mean_KL_values_00 = []
buried_mean_KL_values_003 = []
buried_mean_KL_values_01 = []
buried_mean_KL_values_03 = []
buried_mean_KL_values_06 = []
buried_mean_KL_values_09 = []
buried_mean_KL_values_12 = []

#Mean KL-Divergence for Partially Buried Sites
intermediate_mean_KL_values_00 = []
intermediate_mean_KL_values_003 = []
intermediate_mean_KL_values_01 = []
intermediate_mean_KL_values_03 = []
intermediate_mean_KL_values_06 = []
intermediate_mean_KL_values_09 = []
intermediate_mean_KL_values_12 = []

#Mean KL-Divergence for Surface Sites
surface_mean_KL_values_00 = []
surface_mean_KL_values_003 = []
surface_mean_KL_values_01 = []
surface_mean_KL_values_03 = []
surface_mean_KL_values_06 = []
surface_mean_KL_values_09 = []
surface_mean_KL_values_12 = []

#Mean Entropy for Buried Sites
buried_mean_entropy_values_00 = []
buried_mean_entropy_values_003 = []
buried_mean_entropy_values_01 = []
buried_mean_entropy_values_03 = []
buried_mean_entropy_values_06 = []
buried_mean_entropy_values_09 = []
buried_mean_entropy_values_12 = []

#Mean Entropy for Partially Buried Sites
intermediate_mean_entropy_values_00 = []
intermediate_mean_entropy_values_003 = []
intermediate_mean_entropy_values_01 = []
intermediate_mean_entropy_values_03 = []
intermediate_mean_entropy_values_06 = []
intermediate_mean_entropy_values_09 = []
intermediate_mean_entropy_values_12 = []

#Mean Entropy for Surface Sites
surface_mean_entropy_values_00 = []
surface_mean_entropy_values_003 = []
surface_mean_entropy_values_01 = []
surface_mean_entropy_values_03 = []
surface_mean_entropy_values_06 = []
surface_mean_entropy_values_09 = []
surface_mean_entropy_values_12 = []

buried_natural_mean_KL_values = []
intermediate_natural_mean_KL_values = []
surface_natural_mean_KL_values = []

buried_natural_mean_entropy_values = []
intermediate_natural_mean_entropy_values = []
surface_natural_mean_entropy_values = []


all_entropies_00 = []
all_entropies_003 = []
all_entropies_01 = []
all_entropies_03 = []
all_entropies_06 = []
all_entropies_09 = []
all_entropies_12 = []
all_natural_entropies = []




pdb_names = []
chain_names = []
for path, names, filename in os.walk('.',False): #Searchs for all the *.dat files for the #natural protein 
	for file in filename:
		if(re.search(searchStr, file)!=None):
			fileparts = re.split("_",file)
			pdb_id = fileparts[4].upper() #Gets the PDB Name and the chain_id
			chain_id = fileparts[5]
			chain_id = chain_id[0]
			print "Processsing file: " + file	
			pdb_names.append(pdb_id)
			chain_names.append(chain_id)
			
			natural_proteins = file #Open the files with results
			designed_proteins_00 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.0)  + ".dat"
			designed_proteins_01 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.1)  + ".dat"
			designed_proteins_03 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.3)  + ".dat"
			designed_proteins_06 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.6)  + ".dat"
			designed_proteins_09 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.9)  + ".dat"
			designed_proteins_12 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(1.2)  + ".dat"
			designed_proteins_003 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.03)  + ".dat"
			
			split_natural_1 = "align_natural_sample1_data_array_" + pdb_id + "_" + chain_id + ".dat"
			split_natural_2 = "align_natural_sample2_data_array_" + pdb_id + "_" + chain_id + ".dat"
			
			#Calculates all of the data for comparison (ex. entropy)
			natural_distribution = analysis_functions.get_AA_distribution(natural_proteins)     
			natural_entropy = analysis_functions.get_native_entropy(natural_proteins)
			natural_entropy_array = analysis_functions.make_array(natural_entropy) 
			natural_RSA = analysis_functions.get_RSA_Values(natural_proteins)
			natural_RSA_array = analysis_functions.make_array(natural_RSA)
			natural_mean_RSA_values.append(mean(natural_RSA_array)) 
			natural_mean_entropy_values.append(mean(natural_entropy_array)) 
			
			designed_distribution_00 = analysis_functions.get_AA_distribution(designed_proteins_00)        
			designed_entropy_00 = analysis_functions.get_native_entropy(designed_proteins_00)
			designed_entropy_array_00 = analysis_functions.make_array(designed_entropy_00) 
			designed_RSA_00 = analysis_functions.get_RSA_Values(designed_proteins_00)
			designed_RSA_array_00 = analysis_functions.make_array(designed_RSA_00)
			designed_mean_RSA_values_00.append(mean(designed_RSA_array_00)) 
			designed_mean_entropy_values_00.append(mean(designed_entropy_array_00)) 
			
						
			
			designed_distribution_01 = analysis_functions.get_AA_distribution(designed_proteins_01)
			designed_entropy_01 = analysis_functions.get_native_entropy(designed_proteins_01)
			designed_entropy_array_01 = analysis_functions.make_array(designed_entropy_01) 
			designed_RSA_01 = analysis_functions.get_RSA_Values(designed_proteins_01)
			designed_RSA_array_01 = analysis_functions.make_array(designed_RSA_01)
			designed_mean_RSA_values_01.append(mean(designed_RSA_array_01)) 
			designed_mean_entropy_values_01.append(mean(designed_entropy_array_01)) 
			
			
			designed_distribution_03 = analysis_functions.get_AA_distribution(designed_proteins_03)
			designed_entropy_03 = analysis_functions.get_native_entropy(designed_proteins_03)
			designed_entropy_array_03 = analysis_functions.make_array(designed_entropy_03)
			designed_RSA_03 = analysis_functions.get_RSA_Values(designed_proteins_03)
			designed_RSA_array_03 = analysis_functions.make_array(designed_RSA_03)
			designed_mean_RSA_values_03.append(mean(designed_RSA_array_03)) 
			designed_mean_entropy_values_03.append(mean(designed_entropy_array_03)) 
			
			
			designed_distribution_06 = analysis_functions.get_AA_distribution(designed_proteins_06)
			designed_entropy_06 = analysis_functions.get_native_entropy(designed_proteins_06)
			designed_entropy_array_06 = analysis_functions.make_array(designed_entropy_06)
			designed_RSA_06 = analysis_functions.get_RSA_Values(designed_proteins_06)
			designed_RSA_array_06 = analysis_functions.make_array(designed_RSA_06)
			designed_mean_RSA_values_06.append(mean(designed_RSA_array_06)) 
			designed_mean_entropy_values_06.append(mean(designed_entropy_array_06)) 
			
			
			designed_distribution_09 = analysis_functions.get_AA_distribution(designed_proteins_09)
			designed_entropy_09 = analysis_functions.get_native_entropy(designed_proteins_09)
			designed_entropy_array_09 = analysis_functions.make_array(designed_entropy_09)
			designed_RSA_09 = analysis_functions.get_RSA_Values(designed_proteins_09)
			designed_RSA_array_09 = analysis_functions.make_array(designed_RSA_09)
			designed_mean_RSA_values_09.append(mean(designed_RSA_array_09)) 
			designed_mean_entropy_values_09.append(mean(designed_entropy_array_09)) 
			
			
			designed_distribution_12 = analysis_functions.get_AA_distribution(designed_proteins_12)
			designed_entropy_12 = analysis_functions.get_native_entropy(designed_proteins_12)
			designed_entropy_array_12 = analysis_functions.make_array(designed_entropy_12)
			designed_RSA_12 = analysis_functions.get_RSA_Values(designed_proteins_12)
			designed_RSA_array_12 = analysis_functions.make_array(designed_RSA_12)
			designed_mean_RSA_values_12.append(mean(designed_RSA_array_12)) 
			designed_mean_entropy_values_12.append(mean(designed_entropy_array_12)) 
			
			designed_distribution_003 = analysis_functions.get_AA_distribution(designed_proteins_003)
			designed_entropy_003 = analysis_functions.get_native_entropy(designed_proteins_003)
			designed_entropy_array_003 = analysis_functions.make_array(designed_entropy_003)
			designed_RSA_003 = analysis_functions.get_RSA_Values(designed_proteins_003)
			designed_RSA_array_003 = analysis_functions.make_array(designed_RSA_003)
			designed_mean_RSA_values_003.append(mean(designed_RSA_array_003)) 
			designed_mean_entropy_values_003.append(mean(designed_entropy_array_003)) 
			
			
			[natural_RSA_entropy_corr,p_value] = pearsonr(natural_RSA_array, natural_entropy_array) 
			natural_RSA_entropy_corr = float(natural_RSA_entropy_corr)
			natural_cor_entropy_RSA_values.append(natural_RSA_entropy_corr)
			
			
			[designed_RSA_entropy_corr_003,p_value] = pearsonr(designed_RSA_array_003, designed_entropy_array_003) 
			designed_RSA_entropy_corr_003 = float(designed_RSA_entropy_corr_003)
			designed_cor_entropy_RSA_values_003.append(designed_RSA_entropy_corr_003)
			
			[designed_RSA_entropy_corr_00,p_value] = pearsonr(designed_RSA_array_00, designed_entropy_array_00) 
			designed_RSA_entropy_corr_00 = float(designed_RSA_entropy_corr_00)
			designed_cor_entropy_RSA_values_00.append(designed_RSA_entropy_corr_00)
			
			[designed_RSA_entropy_corr_01,p_value] = pearsonr(designed_RSA_array_01, designed_entropy_array_01) 
			designed_RSA_entropy_corr_01 = float(designed_RSA_entropy_corr_01)
			designed_cor_entropy_RSA_values_01.append(designed_RSA_entropy_corr_01)
			
			[designed_RSA_entropy_corr_03,p_value] = pearsonr(designed_RSA_array_03, designed_entropy_array_03) 
			designed_RSA_entropy_corr_03 = float(designed_RSA_entropy_corr_03)
			designed_cor_entropy_RSA_values_03.append(designed_RSA_entropy_corr_03)
			
			[designed_RSA_entropy_corr_06,p_value] = pearsonr(designed_RSA_array_06, designed_entropy_array_06) 
			designed_RSA_entropy_corr_06 = float(designed_RSA_entropy_corr_06)
			designed_cor_entropy_RSA_values_06.append(designed_RSA_entropy_corr_06)
			
			[designed_RSA_entropy_corr_09,p_value] = pearsonr(designed_RSA_array_09, designed_entropy_array_09) 
			designed_RSA_entropy_corr_09 = float(designed_RSA_entropy_corr_09)
			designed_cor_entropy_RSA_values_09.append(designed_RSA_entropy_corr_09)
			
			[designed_RSA_entropy_corr_12,p_value] = pearsonr(designed_RSA_array_12, designed_entropy_array_12) 
			designed_RSA_entropy_corr_12 = float(designed_RSA_entropy_corr_12)
			designed_cor_entropy_RSA_values_12.append(designed_RSA_entropy_corr_12)
			
			#Calculates the KL Divergence Values 
			KL_00 = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_00)
			KL_01 = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_01)
			KL_03 = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_03)
			KL_06 = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_06)
			KL_09 = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_09)
			KL_12 = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_12)
			KL_003 = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_003)
			
			KL_array_00 = array(KL_00)
			KL_array_003 = array(KL_003)
			KL_array_01 = array(KL_01)
			KL_array_03 = array(KL_03)
			KL_array_06 = array(KL_06)
			KL_array_09 = array(KL_09)
			KL_array_12 = array(KL_12)
			
			designed_mean_KL_values_00.append(mean(KL_array_00))
			designed_mean_KL_values_003.append(mean(KL_array_003))
			designed_mean_KL_values_01.append(mean(KL_array_01))
			designed_mean_KL_values_03.append(mean(KL_array_03))
			designed_mean_KL_values_06.append(mean(KL_array_06))
			designed_mean_KL_values_09.append(mean(KL_array_09))
			designed_mean_KL_values_12.append(mean(KL_array_12))
			
			split_1_distribution = analysis_functions.get_AA_distribution(split_natural_1)
			split_2_distribution = analysis_functions.get_AA_distribution(split_natural_2)
			split_KL = analysis_functions.get_Kullback_Leibler(split_1_distribution, split_2_distribution)
			split_KL_array = array(split_KL)
			natural_mean_split_KL_values.append(mean(split_KL_array))
			'''
			#These lines seperate the data into three classes: buried, partially buried, and surface sites
			[buried_entropy_values_00, buried_KL_values_00, intermediate_entropy_values_00, intermediate_KL_values_00, surface_entropy_values_00, surface_KL_values_00] = analysis_functions.get_position_dependent_data(designed_RSA_00, designed_entropy_00, KL_00)
			
			[buried_entropy_values_003, buried_KL_values_003, intermediate_entropy_values_003, intermediate_KL_values_003, surface_entropy_values_003, surface_KL_values_003] = analysis_functions.get_position_dependent_data(designed_RSA_003, designed_entropy_003, KL_003)
			
			[buried_entropy_values_01, buried_KL_values_01, intermediate_entropy_values_01, intermediate_KL_values_01, surface_entropy_values_01, surface_KL_values_01] = analysis_functions.get_position_dependent_data(designed_RSA_01, designed_entropy_01, KL_01)
			
			[buried_entropy_values_03, buried_KL_values_03, intermediate_entropy_values_03, intermediate_KL_values_03, surface_entropy_values_03, surface_KL_values_03] = analysis_functions.get_position_dependent_data(designed_RSA_03, designed_entropy_03, KL_03)
			
			[buried_entropy_values_06, buried_KL_values_06, intermediate_entropy_values_06, intermediate_KL_values_06, surface_entropy_values_06, surface_KL_values_06] = analysis_functions.get_position_dependent_data(designed_RSA_06, designed_entropy_06, KL_06)
			
			[buried_entropy_values_09, buried_KL_values_09, intermediate_entropy_values_09, intermediate_KL_values_09, surface_entropy_values_09, surface_KL_values_09] = analysis_functions.get_position_dependent_data(designed_RSA_09, designed_entropy_09, KL_09)
			
			[buried_entropy_values_12, buried_KL_values_12, intermediate_entropy_values_12, intermediate_KL_values_12, surface_entropy_values_12, surface_KL_values_12] = analysis_functions.get_position_dependent_data(designed_RSA_12, designed_entropy_12, KL_12)
			
			[natural_buried_entropy_values, natural_buried_KL_values, natural_intermediate_entropy_values, natural_intermediate_KL_values, natural_surface_entropy_values, natural_surface_KL_values] = analysis_functions.get_position_dependent_data(natural_RSA, natural_entropy, split_KL)
			 
			#Store the data that was just calculated 
			buried_KL_array_00 = array(buried_KL_values_00)
			buried_KL_array_003 = array(buried_KL_values_003)
			buried_KL_array_01 = array(buried_KL_values_01)
			buried_KL_array_03 = array(buried_KL_values_03)
			buried_KL_array_06 = array(buried_KL_values_06)
			buried_KL_array_09 = array(buried_KL_values_09)
			buried_KL_array_12 = array(buried_KL_values_12)
			 
			buried_entropy_array_00 = array(buried_entropy_values_00)
			buried_entropy_array_003 = array(buried_entropy_values_003)
			buried_entropy_array_01 = array(buried_entropy_values_01)
			buried_entropy_array_03 = array(buried_entropy_values_03)
			buried_entropy_array_06 = array(buried_entropy_values_06)
			buried_entropy_array_09 = array(buried_entropy_values_09)
			buried_entropy_array_12 = array(buried_entropy_values_12)
			
			buried_mean_entropy_values_00.append(mean(buried_entropy_array_00))
			buried_mean_entropy_values_003.append(mean(buried_entropy_array_003))
			buried_mean_entropy_values_01.append(mean(buried_entropy_array_01))
			buried_mean_entropy_values_03.append(mean(buried_entropy_array_03))
			buried_mean_entropy_values_06.append(mean(buried_entropy_array_06))
			buried_mean_entropy_values_09.append(mean(buried_entropy_array_09))
			buried_mean_entropy_values_12.append(mean(buried_entropy_array_12))  
			
			buried_mean_KL_values_00.append(mean(array(buried_KL_array_00)))
			buried_mean_KL_values_003.append(mean(array(buried_KL_array_003)))
			buried_mean_KL_values_01.append(mean(array(buried_KL_array_01)))
			buried_mean_KL_values_03.append(mean(array(buried_KL_array_03)))
			buried_mean_KL_values_06.append(mean(array(buried_KL_array_06)))
			buried_mean_KL_values_09.append(mean(array(buried_KL_array_09)))
			buried_mean_KL_values_12.append(mean(array(buried_KL_array_12)))
			
			intermediate_KL_array_00 = array(intermediate_KL_values_00)
			intermediate_KL_array_003 = array(intermediate_KL_values_003)
			intermediate_KL_array_01 = array(intermediate_KL_values_01)
			intermediate_KL_array_03 = array(intermediate_KL_values_03)
			intermediate_KL_array_06 = array(intermediate_KL_values_06)
			intermediate_KL_array_09 = array(intermediate_KL_values_09)
			intermediate_KL_array_12 = array(intermediate_KL_values_12)
			
			intermediate_entropy_array_00 = array(intermediate_entropy_values_00)
			intermediate_entropy_array_003 = array(intermediate_entropy_values_003)
			intermediate_entropy_array_01 = array(intermediate_entropy_values_01)
			intermediate_entropy_array_03 = array(intermediate_entropy_values_03)
			intermediate_entropy_array_06 = array(intermediate_entropy_values_06)
			intermediate_entropy_array_09 = array(intermediate_entropy_values_09)
			intermediate_entropy_array_12 = array(intermediate_entropy_values_12)
			
			intermediate_mean_entropy_values_00.append(mean(intermediate_entropy_array_00))
			intermediate_mean_entropy_values_003.append(mean(intermediate_entropy_array_003))
			intermediate_mean_entropy_values_01.append(mean(intermediate_entropy_array_01))
			intermediate_mean_entropy_values_03.append(mean(intermediate_entropy_array_03))
			intermediate_mean_entropy_values_06.append(mean(intermediate_entropy_array_06))
			intermediate_mean_entropy_values_09.append(mean(intermediate_entropy_array_09))
			intermediate_mean_entropy_values_12.append(mean(intermediate_entropy_array_12))
			
			intermediate_mean_KL_values_00.append(mean(array(intermediate_KL_array_00)))  
			intermediate_mean_KL_values_003.append(mean(array(intermediate_KL_array_003))) 
			intermediate_mean_KL_values_01.append(mean(array(intermediate_KL_array_01))) 
			intermediate_mean_KL_values_03.append(mean(array(intermediate_KL_array_03))) 
			intermediate_mean_KL_values_06.append(mean(array(intermediate_KL_array_06))) 
			intermediate_mean_KL_values_09.append(mean(array(intermediate_KL_array_09))) 
			intermediate_mean_KL_values_12.append(mean(array(intermediate_KL_array_12))) 
				
			surface_KL_array_00 = array(surface_KL_values_00)
			surface_KL_array_003 = array(surface_KL_values_003)
			surface_KL_array_01 = array(surface_KL_values_01)
			surface_KL_array_03 = array(surface_KL_values_03)
			surface_KL_array_06 = array(surface_KL_values_06)
			surface_KL_array_09 = array(surface_KL_values_09)
			surface_KL_array_12 = array(surface_KL_values_12)
			
			surface_entropy_array_00 = array(surface_entropy_values_00)
			surface_entropy_array_003 = array(surface_entropy_values_003)
			surface_entropy_array_01 = array(surface_entropy_values_01)
			surface_entropy_array_03 = array(surface_entropy_values_03)
			surface_entropy_array_06 = array(surface_entropy_values_06)
			surface_entropy_array_09 = array(surface_entropy_values_09)
			surface_entropy_array_12 = array(surface_entropy_values_12)
			
			surface_mean_KL_values_00.append(mean(array(surface_KL_array_00)))
			surface_mean_KL_values_003.append(mean(array(surface_KL_array_003)))
			surface_mean_KL_values_01.append(mean(array(surface_KL_array_01)))
			surface_mean_KL_values_03.append(mean(array(surface_KL_array_03)))
			surface_mean_KL_values_06.append(mean(array(surface_KL_array_06)))
			surface_mean_KL_values_09.append(mean(array(surface_KL_array_09)))
			surface_mean_KL_values_12.append(mean(array(surface_KL_array_12)))
			
			surface_mean_entropy_values_00.append(mean(surface_entropy_array_00))
			surface_mean_entropy_values_003.append(mean(surface_entropy_array_003))
			surface_mean_entropy_values_01.append(mean(surface_entropy_array_01))
			surface_mean_entropy_values_03.append(mean(surface_entropy_array_03))
			surface_mean_entropy_values_06.append(mean(surface_entropy_array_06))
			surface_mean_entropy_values_09.append(mean(surface_entropy_array_09))
			surface_mean_entropy_values_12.append(mean(surface_entropy_array_12))
			
			buried_natural_mean_KL_values.append(mean(array(natural_buried_KL_values)))
			intermediate_natural_mean_KL_values.append(mean(array(natural_intermediate_KL_values)))
			surface_natural_mean_KL_values.append(mean(array(natural_surface_KL_values)))
			
			buried_natural_mean_entropy_values.append(mean(array(natural_buried_entropy_values)))
			intermediate_natural_mean_entropy_values.append(mean(array(natural_intermediate_entropy_values)))
			surface_natural_mean_entropy_values.append(mean(array(natural_surface_entropy_values)))			
			'''
			
			for element in designed_entropy_00:
				if element!= '\n':
					all_entropies_00.append(element)		
			#print all_entropies_00
			#print "\n"

			for element in designed_entropy_06:
				if element!= '\n':
					all_entropies_06.append(element)		
			#print all_entropies_12
			#print "\n"
			
			for element in natural_entropy:
				if element!= '\n':
					all_natural_entropies.append(element)		
			#print all_natural_entropies
						

all_entropies_00_array = analysis_functions.make_array(all_entropies_00)
all_entropies_06_array = analysis_functions.make_array(all_entropies_06)
all_natural_entropies_array = analysis_functions.make_array(all_natural_entropies)

[fixed_natural_t, fixed_natural_pvalue] = ttest_rel(all_entropies_00_array, all_natural_entropies_array)

print "T-stat:  ", fixed_natural_t
print "p-value: ", fixed_natural_pvalue

			
			
			
			
			
			
			
			
			
			
			
			