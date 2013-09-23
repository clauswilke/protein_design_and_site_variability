import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from scipy.stats import pearsonr
import analysis_functions as af

#Last Date Updated: June 24, 2013
#Description: This is a script that takes all of the *.dat files with the generated AA count and RSA values and then calculates the necessary values for comparision. Ex. Entropy, KL Divergence

#This is dictionary that maps the natural sequence identifier with its PDB for the Noah Dataset
identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076':'2X1B', \
	'PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37', \
	'PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY', 'PF00355' :'1FQT', \
	'PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS', 'PF00542' :'1CTF', \
	'PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8',  \
	'PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI', 'PF02823' :'1AQT', \
	'PF04002':'2QLC','PF07686' :'2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

#Get files 
#This searches for all of the unaligned sequences 
PDBS = identity_dict.values()
#chStr = "align_natural_data_array_" + "[a-zA-Z0-9_\.\-]*" +  ".dat"
#find all csv files that match the search string
temp = 0.0
temp_array = [0.0, 0.1, 0.3, 0.6, 0.9, 1.2]
data = []
count = 1 
natural_mean_RSA_values = []
natural_mean_entropy_values = []
natural_cor_entropy_RSA_values = []

designed_mean_RSA_values_00 = []
designed_mean_RSA_values_03 = []
designed_mean_RSA_values_06 = []
designed_mean_RSA_values_09 = []
designed_mean_RSA_values_12 = []
designed_mean_RSA_values_18 = []
designed_mean_RSA_values_24 = []
designed_mean_RSA_values_soft = []

natural_mean_RSA_values = []
natural_mean_entropy_values = []

KL_array_00 = []
KL_array_03 = []
KL_array_06 = []
KL_array_09 = []
KL_array_12 = []
KL_array_18 = []
KL_array_24 = []
KL_array_soft = []

designed_mean_KL_values_00 = []
designed_mean_KL_values_03 = []
designed_mean_KL_values_06 = []
designed_mean_KL_values_09 = []
designed_mean_KL_values_12 = []
designed_mean_KL_values_18 = []
designed_mean_KL_values_24 = []
designed_mean_KL_values_soft = []

designed_mean_entropy_values_00 = []
designed_mean_entropy_values_03 = []
designed_mean_entropy_values_06 = []
designed_mean_entropy_values_09 = []
designed_mean_entropy_values_12 = []
designed_mean_entropy_values_18 = []
designed_mean_entropy_values_24 = []
designed_mean_entropy_values_soft = []

designed_cor_entropy_RSA_values_00 = []
designed_cor_entropy_RSA_values_03 = []
designed_cor_entropy_RSA_values_06 = []
designed_cor_entropy_RSA_values_09 = []
designed_cor_entropy_RSA_values_12 = []
designed_cor_entropy_RSA_values_18 = []
designed_cor_entropy_RSA_values_24 = []
designed_cor_entropy_RSA_values_soft = []

split_KL_array = []
natural_mean_split_KL_values = []

buried_mean_KL_values_00 = []
buried_mean_KL_values_03 = []
buried_mean_KL_values_06 = []
buried_mean_KL_values_09 = []
buried_mean_KL_values_12 = []
buried_mean_KL_values_18 = []
buried_mean_KL_values_24 = []
buried_mean_KL_values_soft = []

intermediate_mean_KL_values_00 = []
intermediate_mean_KL_values_03 = []
intermediate_mean_KL_values_06 = []
intermediate_mean_KL_values_09 = []
intermediate_mean_KL_values_12 = []
intermediate_mean_KL_values_18 = []
intermediate_mean_KL_values_24 = []
intermediate_mean_KL_values_soft = []

surface_mean_KL_values_00 = []
surface_mean_KL_values_03 = []
surface_mean_KL_values_06 = []
surface_mean_KL_values_09 = []
surface_mean_KL_values_12 = []
surface_mean_KL_values_18 = []
surface_mean_KL_values_24 = []
surface_mean_KL_values_soft = []

buried_mean_entropy_values_00 = []
buried_mean_entropy_values_03 = []
buried_mean_entropy_values_06 = []
buried_mean_entropy_values_09 = []
buried_mean_entropy_values_12 = []
buried_mean_entropy_values_18 = []
buried_mean_entropy_values_24 = []
buried_mean_entropy_values_soft = []

intermediate_mean_entropy_values_00 = []
intermediate_mean_entropy_values_03 = []
intermediate_mean_entropy_values_06 = []
intermediate_mean_entropy_values_09 = []
intermediate_mean_entropy_values_12 = []
intermediate_mean_entropy_values_18 = []
intermediate_mean_entropy_values_24 = []
intermediate_mean_entropy_values_soft = []

surface_mean_entropy_values_00 = []
surface_mean_entropy_values_03 = []
surface_mean_entropy_values_06 = []
surface_mean_entropy_values_09 = []
surface_mean_entropy_values_12 = []
surface_mean_entropy_values_18 = []
surface_mean_entropy_values_24 = []
surface_mean_entropy_values_soft = []

buried_natural_mean_KL_values = []
intermediate_natural_mean_KL_values = []
surface_natural_mean_KL_values = []

buried_natural_mean_entropy_values = []
intermediate_natural_mean_entropy_values = []
surface_natural_mean_entropy_values = []

pdb_names = []
chain_names = []
chain_id = "A"
for pdb_id in PDBS:
	file = "align_natural_data_array_" + pdb_id + "_" + chain_id +  ".dat"	
	print "Processsing file: " + file	
	pdb_names.append(pdb_id)
	chain_names.append(chain_id)
	
	natural_proteins = file
	designed_proteins_00 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.0)  + ".dat"
	designed_proteins_03 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.3)  + ".dat"
	designed_proteins_06 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.6)  + ".dat"
	designed_proteins_09 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(0.9)  + ".dat"
	designed_proteins_12 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(1.2)  + ".dat"
	designed_proteins_18 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(1.8)  + ".dat"
	designed_proteins_24 = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str(2.4)  + ".dat"
	designed_proteins_soft = "align_data_array_" + pdb_id + "_" + chain_id + "_" + "soft"  + ".dat"
	
	split_natural_1 = "align_natural_sample1_data_array_" + pdb_id + "_" + chain_id + ".dat"
	split_natural_2 = "align_natural_sample2_data_array_" + pdb_id + "_" + chain_id + ".dat"
	
	#Calculates all of the data for comparison (ex. entropy)
	natural_distribution = af.get_AA_distribution_mod(natural_proteins)     
	natural_entropy = af.get_native_entropy(natural_proteins)
	natural_entropy_array = af.make_array(natural_entropy) 
	natural_RSA = af.get_RSA_Values(natural_proteins)
	natural_RSA_array = af.make_array(natural_RSA)
	natural_mean_RSA_values.append(mean(natural_RSA_array)) 
	natural_mean_entropy_values.append(mean(natural_entropy_array)) 
	
	designed_distribution_00 = af.get_AA_distribution_mod(designed_proteins_00)        
	designed_entropy_00 = af.get_native_entropy(designed_proteins_00)
	designed_entropy_array_00 = af.make_array(designed_entropy_00) 
	designed_RSA_00 = af.get_RSA_Values(designed_proteins_00)
	designed_RSA_array_00 = af.make_array(designed_RSA_00)
	designed_mean_RSA_values_00.append(mean(designed_RSA_array_00)) 
	designed_mean_entropy_values_00.append(mean(designed_entropy_array_00)) 
	
	designed_distribution_03 = af.get_AA_distribution_mod(designed_proteins_03)
	designed_entropy_03 = af.get_native_entropy(designed_proteins_03)
	designed_entropy_array_03 = af.make_array(designed_entropy_03)
	designed_RSA_03 = af.get_RSA_Values(designed_proteins_03)
	designed_RSA_array_03 = af.make_array(designed_RSA_03)
	designed_mean_RSA_values_03.append(mean(designed_RSA_array_03)) 
	designed_mean_entropy_values_03.append(mean(designed_entropy_array_03)) 
	
	designed_distribution_06 = af.get_AA_distribution_mod(designed_proteins_06)
	designed_entropy_06 = af.get_native_entropy(designed_proteins_06)
	designed_entropy_array_06 = af.make_array(designed_entropy_06)
	designed_RSA_06 = af.get_RSA_Values(designed_proteins_06)
	designed_RSA_array_06 = af.make_array(designed_RSA_06)
	designed_mean_RSA_values_06.append(mean(designed_RSA_array_06)) 
	designed_mean_entropy_values_06.append(mean(designed_entropy_array_06)) 
	
	designed_distribution_09 = af.get_AA_distribution_mod(designed_proteins_09)
	designed_entropy_09 = af.get_native_entropy(designed_proteins_09)
	designed_entropy_array_09 = af.make_array(designed_entropy_09)
	designed_RSA_09 = af.get_RSA_Values(designed_proteins_09)
	designed_RSA_array_09 = af.make_array(designed_RSA_09)
	designed_mean_RSA_values_09.append(mean(designed_RSA_array_09)) 
	designed_mean_entropy_values_09.append(mean(designed_entropy_array_09)) 
	
	designed_distribution_12 = af.get_AA_distribution_mod(designed_proteins_12)
	designed_entropy_12 = af.get_native_entropy(designed_proteins_12)
	designed_entropy_array_12 = af.make_array(designed_entropy_12)
	designed_RSA_12 = af.get_RSA_Values(designed_proteins_12)
	designed_RSA_array_12 = af.make_array(designed_RSA_12)
	designed_mean_RSA_values_12.append(mean(designed_RSA_array_12)) 
	designed_mean_entropy_values_12.append(mean(designed_entropy_array_12)) 
	
	designed_distribution_18 = af.get_AA_distribution_mod(designed_proteins_18)
	designed_entropy_18 = af.get_native_entropy(designed_proteins_18)
	designed_entropy_array_18 = af.make_array(designed_entropy_18)
	designed_RSA_18 = af.get_RSA_Values(designed_proteins_18)
	designed_RSA_array_18 = af.make_array(designed_RSA_18)
	designed_mean_RSA_values_18.append(mean(designed_RSA_array_18)) 
	designed_mean_entropy_values_18.append(mean(designed_entropy_array_18)) 
	
	designed_distribution_24 = af.get_AA_distribution_mod(designed_proteins_24)
	designed_entropy_24 = af.get_native_entropy(designed_proteins_24)
	designed_entropy_array_24 = af.make_array(designed_entropy_24)
	designed_RSA_24 = af.get_RSA_Values(designed_proteins_24)
	designed_RSA_array_24 = af.make_array(designed_RSA_24)
	designed_mean_RSA_values_24.append(mean(designed_RSA_array_24)) 
	designed_mean_entropy_values_24.append(mean(designed_entropy_array_24)) 
	
	designed_distribution_soft = af.get_AA_distribution_mod(designed_proteins_soft)
	designed_entropy_soft = af.get_native_entropy(designed_proteins_soft)
	designed_entropy_array_soft = af.make_array(designed_entropy_soft)
	designed_RSA_soft = af.get_RSA_Values(designed_proteins_soft)
	designed_RSA_array_soft = af.make_array(designed_RSA_soft)
	designed_mean_RSA_values_soft.append(mean(designed_RSA_array_soft)) 
	designed_mean_entropy_values_soft.append(mean(designed_entropy_array_soft)) 
	
	[natural_RSA_entropy_corr,p_value] = pearsonr(natural_RSA_array, natural_entropy_array) 
	natural_RSA_entropy_corr = float(natural_RSA_entropy_corr)
	natural_cor_entropy_RSA_values.append(natural_RSA_entropy_corr)
	
	[designed_RSA_entropy_corr_00,p_value] = pearsonr(designed_RSA_array_00, designed_entropy_array_00) 
	designed_RSA_entropy_corr_00 = float(designed_RSA_entropy_corr_00)
	designed_cor_entropy_RSA_values_00.append(designed_RSA_entropy_corr_00)
	
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
	
	[designed_RSA_entropy_corr_18,p_value] = pearsonr(designed_RSA_array_18, designed_entropy_array_18) 
	designed_RSA_entropy_corr_18 = float(designed_RSA_entropy_corr_18)
	designed_cor_entropy_RSA_values_18.append(designed_RSA_entropy_corr_18)
	
	[designed_RSA_entropy_corr_18,p_value] = pearsonr(designed_RSA_array_18, designed_entropy_array_18) 
	designed_RSA_entropy_corr_18 = float(designed_RSA_entropy_corr_18)
	designed_cor_entropy_RSA_values_18.append(designed_RSA_entropy_corr_18)
	
	[designed_RSA_entropy_corr_18,p_value] = pearsonr(designed_RSA_array_18, designed_entropy_array_18) 
	designed_RSA_entropy_corr_18 = float(designed_RSA_entropy_corr_18)
	designed_cor_entropy_RSA_values_18.append(designed_RSA_entropy_corr_18)
	
	[designed_RSA_entropy_corr_24,p_value] = pearsonr(designed_RSA_array_24, designed_entropy_array_24) 
	designed_RSA_entropy_corr_24 = float(designed_RSA_entropy_corr_24)
	designed_cor_entropy_RSA_values_24.append(designed_RSA_entropy_corr_24)
	
	[designed_RSA_entropy_corr_soft,p_value] = pearsonr(designed_RSA_array_soft, designed_entropy_array_soft) 
	designed_RSA_entropy_corr_soft = float(designed_RSA_entropy_corr_soft)
	designed_cor_entropy_RSA_values_soft.append(designed_RSA_entropy_corr_soft)
	
	#Calculates the KL Divergence Values 
	KL_00 = af.get_Kullback_Leibler(natural_distribution, designed_distribution_00)
	KL_03 = af.get_Kullback_Leibler(natural_distribution, designed_distribution_03)
	KL_06 = af.get_Kullback_Leibler(natural_distribution, designed_distribution_06)
	KL_09 = af.get_Kullback_Leibler(natural_distribution, designed_distribution_09)
	KL_12 = af.get_Kullback_Leibler(natural_distribution, designed_distribution_12)
	KL_18 = af.get_Kullback_Leibler(natural_distribution, designed_distribution_18)
	KL_24 = af.get_Kullback_Leibler(natural_distribution, designed_distribution_24)
	KL_soft = af.get_Kullback_Leibler(natural_distribution, designed_distribution_soft)
	
	KL_array_00 = array(KL_00)
	KL_array_03 = array(KL_03)
	KL_array_06 = array(KL_06)
	KL_array_09 = array(KL_09)
	KL_array_12 = array(KL_12)
	KL_array_18 = array(KL_18)
	KL_array_24 = array(KL_24)
	KL_array_soft = array(KL_soft)
	
	designed_mean_KL_values_00.append(mean(KL_array_00))
	designed_mean_KL_values_03.append(mean(KL_array_03))
	designed_mean_KL_values_06.append(mean(KL_array_06))
	designed_mean_KL_values_09.append(mean(KL_array_09))
	designed_mean_KL_values_12.append(mean(KL_array_12))
	designed_mean_KL_values_18.append(mean(KL_array_18))
	designed_mean_KL_values_24.append(mean(KL_array_24))
	designed_mean_KL_values_soft.append(mean(KL_array_soft))
	
	split_1_distribution = af.get_AA_distribution(split_natural_1)
	split_2_distribution = af.get_AA_distribution(split_natural_2)
	split_KL = af.get_Kullback_Leibler(split_1_distribution, split_2_distribution)
	split_KL_array = array(split_KL)
	natural_mean_split_KL_values.append(mean(split_KL_array))
	
	#These lines seperate the data into three classes: buried, partially buried, and surface sites
	[buried_entropy_values_00, buried_KL_values_00, intermediate_entropy_values_00, intermediate_KL_values_00, surface_entropy_values_00, surface_KL_values_00] = af.get_position_dependent_data(designed_RSA_00, designed_entropy_00, KL_00)
	
	[buried_entropy_values_03, buried_KL_values_03, intermediate_entropy_values_03, intermediate_KL_values_03, surface_entropy_values_03, surface_KL_values_03] = af.get_position_dependent_data(designed_RSA_03, designed_entropy_03, KL_03)
	
	[buried_entropy_values_06, buried_KL_values_06, intermediate_entropy_values_06, intermediate_KL_values_06, surface_entropy_values_06, surface_KL_values_06] = af.get_position_dependent_data(designed_RSA_06, designed_entropy_06, KL_06)
	
	[buried_entropy_values_09, buried_KL_values_09, intermediate_entropy_values_09, intermediate_KL_values_09, surface_entropy_values_09, surface_KL_values_09] = af.get_position_dependent_data(designed_RSA_09, designed_entropy_09, KL_09)
	
	[buried_entropy_values_12, buried_KL_values_12, intermediate_entropy_values_12, intermediate_KL_values_12, surface_entropy_values_12, surface_KL_values_12] = af.get_position_dependent_data(designed_RSA_12, designed_entropy_12, KL_12)
	
	[buried_entropy_values_18, buried_KL_values_18, intermediate_entropy_values_18, intermediate_KL_values_18, surface_entropy_values_18, surface_KL_values_18] = af.get_position_dependent_data(designed_RSA_18, designed_entropy_18, KL_18)
	
	[buried_entropy_values_24, buried_KL_values_24, intermediate_entropy_values_24, intermediate_KL_values_24, surface_entropy_values_24, surface_KL_values_24] = af.get_position_dependent_data(designed_RSA_24, designed_entropy_24, KL_24)
	
	[buried_entropy_values_soft, buried_KL_values_soft, intermediate_entropy_values_soft, intermediate_KL_values_soft, surface_entropy_values_soft, surface_KL_values_soft] = af.get_position_dependent_data(designed_RSA_soft, designed_entropy_soft, KL_soft)
	
	[natural_buried_entropy_values, natural_buried_KL_values, natural_intermediate_entropy_values, natural_intermediate_KL_values, natural_surface_entropy_values, natural_surface_KL_values] = af.get_position_dependent_data(natural_RSA, natural_entropy, split_KL)
	 
	buried_KL_array_00 = array(buried_KL_values_00)
	buried_KL_array_03 = array(buried_KL_values_03)
	buried_KL_array_06 = array(buried_KL_values_06)
	buried_KL_array_09 = array(buried_KL_values_09)
	buried_KL_array_12 = array(buried_KL_values_12)
	buried_KL_array_18 = array(buried_KL_values_18)
	buried_KL_array_24 = array(buried_KL_values_24)
	buried_KL_array_soft = array(buried_KL_values_soft)             
	
	buried_entropy_array_00 = array(buried_entropy_values_00)
	buried_entropy_array_03 = array(buried_entropy_values_03)
	buried_entropy_array_06 = array(buried_entropy_values_06)
	buried_entropy_array_09 = array(buried_entropy_values_09)
	buried_entropy_array_12 = array(buried_entropy_values_12)
	buried_entropy_array_18 = array(buried_entropy_values_18)
	buried_entropy_array_24 = array(buried_entropy_values_24)
	buried_entropy_array_soft = array(buried_entropy_values_soft)
	
	buried_mean_entropy_values_00.append(mean(buried_entropy_array_00))
	buried_mean_entropy_values_03.append(mean(buried_entropy_array_03))
	buried_mean_entropy_values_06.append(mean(buried_entropy_array_06))
	buried_mean_entropy_values_09.append(mean(buried_entropy_array_09))
	buried_mean_entropy_values_12.append(mean(buried_entropy_array_12))  
	buried_mean_entropy_values_18.append(mean(buried_entropy_array_18))
	buried_mean_entropy_values_24.append(mean(buried_entropy_array_24))
	buried_mean_entropy_values_soft.append(mean(buried_entropy_array_soft))
	
	buried_mean_KL_values_00.append(mean(array(buried_KL_array_00)))
	buried_mean_KL_values_03.append(mean(array(buried_KL_array_03)))
	buried_mean_KL_values_06.append(mean(array(buried_KL_array_06)))
	buried_mean_KL_values_09.append(mean(array(buried_KL_array_09)))
	buried_mean_KL_values_12.append(mean(array(buried_KL_array_12)))
	buried_mean_KL_values_18.append(mean(array(buried_KL_array_18)))
	buried_mean_KL_values_24.append(mean(array(buried_KL_array_24)))
	buried_mean_KL_values_soft.append(mean(array(buried_KL_array_soft)))
	
	intermediate_KL_array_00 = array(intermediate_KL_values_00)
	intermediate_KL_array_03 = array(intermediate_KL_values_03)
	intermediate_KL_array_06 = array(intermediate_KL_values_06)
	intermediate_KL_array_09 = array(intermediate_KL_values_09)
	intermediate_KL_array_12 = array(intermediate_KL_values_12)
	intermediate_KL_array_18 = array(intermediate_KL_values_18)
	intermediate_KL_array_24 = array(intermediate_KL_values_24)
	intermediate_KL_array_soft = array(intermediate_KL_values_soft)
	
	intermediate_entropy_array_00 = array(intermediate_entropy_values_00)
	intermediate_entropy_array_03 = array(intermediate_entropy_values_03)
	intermediate_entropy_array_06 = array(intermediate_entropy_values_06)
	intermediate_entropy_array_09 = array(intermediate_entropy_values_09)
	intermediate_entropy_array_12 = array(intermediate_entropy_values_12)
	intermediate_entropy_array_18 = array(intermediate_entropy_values_18)
	intermediate_entropy_array_24 = array(intermediate_entropy_values_24)
	intermediate_entropy_array_soft = array(intermediate_entropy_values_soft)
	
	intermediate_mean_entropy_values_00.append(mean(intermediate_entropy_array_00))
	intermediate_mean_entropy_values_03.append(mean(intermediate_entropy_array_03))
	intermediate_mean_entropy_values_06.append(mean(intermediate_entropy_array_06))
	intermediate_mean_entropy_values_09.append(mean(intermediate_entropy_array_09))
	intermediate_mean_entropy_values_12.append(mean(intermediate_entropy_array_12))
	intermediate_mean_entropy_values_18.append(mean(intermediate_entropy_array_18))           
	intermediate_mean_entropy_values_24.append(mean(intermediate_entropy_array_24))
	intermediate_mean_entropy_values_soft.append(mean(intermediate_entropy_array_soft))
	
	intermediate_mean_KL_values_00.append(mean(array(intermediate_KL_array_00)))  
	intermediate_mean_KL_values_03.append(mean(array(intermediate_KL_array_03))) 
	intermediate_mean_KL_values_06.append(mean(array(intermediate_KL_array_06))) 
	intermediate_mean_KL_values_09.append(mean(array(intermediate_KL_array_09))) 
	intermediate_mean_KL_values_12.append(mean(array(intermediate_KL_array_12))) 
	intermediate_mean_KL_values_18.append(mean(array(intermediate_KL_array_18))) 
	intermediate_mean_KL_values_24.append(mean(array(intermediate_KL_array_24))) 
	intermediate_mean_KL_values_soft.append(mean(array(intermediate_KL_array_soft))) 
	
	surface_KL_array_00 = array(surface_KL_values_00)
	surface_KL_array_03 = array(surface_KL_values_03)
	surface_KL_array_06 = array(surface_KL_values_06)
	surface_KL_array_09 = array(surface_KL_values_09)
	surface_KL_array_12 = array(surface_KL_values_12)
	surface_KL_array_18 = array(surface_KL_values_18)
	surface_KL_array_24 = array(surface_KL_values_24)
	surface_KL_array_soft = array(surface_KL_values_soft)
	
	surface_entropy_array_00 = array(surface_entropy_values_00)
	surface_entropy_array_03 = array(surface_entropy_values_03)
	surface_entropy_array_06 = array(surface_entropy_values_06)
	surface_entropy_array_09 = array(surface_entropy_values_09)
	surface_entropy_array_12 = array(surface_entropy_values_12)
	surface_entropy_array_18 = array(surface_entropy_values_18)
	surface_entropy_array_24 = array(surface_entropy_values_24)
	surface_entropy_array_soft = array(surface_entropy_values_soft)
	
	surface_mean_KL_values_00.append(mean(array(surface_KL_array_00)))
	surface_mean_KL_values_03.append(mean(array(surface_KL_array_03)))
	surface_mean_KL_values_06.append(mean(array(surface_KL_array_06)))
	surface_mean_KL_values_09.append(mean(array(surface_KL_array_09)))
	surface_mean_KL_values_12.append(mean(array(surface_KL_array_12)))
	surface_mean_KL_values_18.append(mean(array(surface_KL_array_18)))
	surface_mean_KL_values_24.append(mean(array(surface_KL_array_24)))
	surface_mean_KL_values_soft.append(mean(array(surface_KL_array_soft)))
	
	surface_mean_entropy_values_00.append(mean(surface_entropy_array_00))
	surface_mean_entropy_values_03.append(mean(surface_entropy_array_03))
	surface_mean_entropy_values_06.append(mean(surface_entropy_array_06))
	surface_mean_entropy_values_09.append(mean(surface_entropy_array_09))
	surface_mean_entropy_values_12.append(mean(surface_entropy_array_12))
	surface_mean_entropy_values_18.append(mean(surface_entropy_array_18))
	surface_mean_entropy_values_24.append(mean(surface_entropy_array_24))
	surface_mean_entropy_values_soft.append(mean(surface_entropy_array_soft))
	
	buried_natural_mean_KL_values.append(mean(array(natural_buried_KL_values)))
	intermediate_natural_mean_KL_values.append(mean(array(natural_intermediate_KL_values)))
	surface_natural_mean_KL_values.append(mean(array(natural_surface_KL_values)))
	
	buried_natural_mean_entropy_values.append(mean(array(natural_buried_entropy_values)))
	intermediate_natural_mean_entropy_values.append(mean(array(natural_intermediate_entropy_values)))
	surface_natural_mean_entropy_values.append(mean(array(natural_surface_entropy_values)))
	
	#open file name
	natural_graphing_filename = "graph_data_" + pdb_id + "_" + chain_id + "_natural_noah" + ".csv"
	natural_file = open(natural_graphing_filename, "w")
	natural_file.write("entropy\t" + af.dump_csv_line3(natural_entropy_array))
	natural_file.write("RSA\t" + af.dump_csv_line3(natural_RSA_array))
	natural_file.close()

	temp_array = [0.0, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4] 
	for temp in temp_array: #These lines write the calculated site data to to the graph data files for each temperature treatment
		graphing_filename = "graph_data_" + pdb_id + "_" + chain_id + "_" + str(temp) + "_noah.csv"
		graph_file = open(graphing_filename, "w")
		if temp == 0.0:      
			graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_00))
			graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_00))
			graph_file.write("KL\t" + af.dump_csv_line3(KL_array_00))

		elif temp == 0.3: #
			graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_03))
			graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_03))
			graph_file.write("KL\t" + af.dump_csv_line3(KL_array_03))
	 
		elif temp == 0.6: #
			graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_06))
			graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_06))
			graph_file.write("KL\t" + af.dump_csv_line3(KL_array_06))

		elif temp == 0.9: # 
			graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_09))
			graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_09))
			graph_file.write("KL\t" + af.dump_csv_line3(KL_array_00))
  
		elif temp == 1.2: #
			graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_12))
			graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_12))
			graph_file.write("KL\t" + af.dump_csv_line3(KL_array_12))

		elif temp == 1.8: #
			graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_18))
			graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_18))
			graph_file.write("KL\t" + af.dump_csv_line3(KL_array_18))

		elif temp == 2.4: #
			graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_24))
			graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_24))
			graph_file.write("KL\t" + af.dump_csv_line3(KL_array_24))

		graph_file.close()

	graphing_filename = "graph_data_" + pdb_id + "_" + chain_id + "_" + "soft_noah" + ".csv"
	graph_file = open(graphing_filename, "w")
	graph_file.write("entropy\t" + af.dump_csv_line3(designed_entropy_array_soft))
	graph_file.write("RSA\t" + af.dump_csv_line3(designed_RSA_array_soft))
	graph_file.write("KL\t" + af.dump_csv_line3(KL_array_soft))

	graph_file.close()


#These files print the MEAN DATA FOR EACH PDB (Seperated by tabs)
natural_mean_RSA_values_array = array(natural_mean_RSA_values)
natural_mean_entropy_values_array = array(natural_mean_entropy_values)
natural_cor_entropy_RSA_values_array = array(natural_cor_entropy_RSA_values)       
natural_mean_split_KL_array = array(natural_mean_split_KL_values)

natural_mean_buried_KL_array = array(buried_natural_mean_KL_values)
natural_mean_intermediate_KL_array = array(intermediate_natural_mean_KL_values)
natural_mean_surface_KL_array = array(surface_natural_mean_KL_values)

natural_mean_buried_entropy_array = array(buried_natural_mean_entropy_values)
natural_mean_intermediate_entropy_array = array(intermediate_natural_mean_entropy_values)
natural_mean_surface_entropy_array = array(surface_natural_mean_entropy_values)

#open file name
natural_mean_graphing_filename = "graph_mean_data_natural_noah.csv"
natural_mean_file = open(natural_mean_graphing_filename, "w")
natural_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tsplit_mean_KL\n")
natural_name_length = len(pdb_names)
length_counter = 0
while(length_counter<natural_name_length):
    natural_mean_RSA = natural_mean_RSA_values[length_counter]
    natural_mean_entropy = natural_mean_entropy_values[length_counter]
    natural_cor_entropy_RSA = natural_cor_entropy_RSA_values[length_counter] 
    natural_mean_split_KL = natural_mean_split_KL_array[length_counter]
    natural_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(natural_mean_RSA) + "\t" + str(natural_mean_entropy) + "\t" + str(natural_cor_entropy_RSA) + "\t" + str(natural_mean_split_KL) + "\n"
    natural_mean_file.write(natural_filestring)
    length_counter = length_counter + 1
natural_mean_file.close()

#open file name
natural_mean_graphing_filename = "graph_mean_buried_data_natural_noah.csv"
natural_mean_file = open(natural_mean_graphing_filename, "w")
natural_mean_file.write("PDB\tchain\tmean_entropy\tsplit_mean_KL\n")
natural_name_length = len(pdb_names)
length_counter = 0
while(length_counter<natural_name_length):
    buried_natural_mean_entropy = natural_mean_buried_entropy_array[length_counter]
    buried_natural_mean_split_KL = natural_mean_buried_KL_array[length_counter]
    natural_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(buried_natural_mean_entropy) + "\t" + str(natural_mean_split_KL) + "\n"
    natural_mean_file.write(natural_filestring)
    length_counter = length_counter + 1
natural_mean_file.close()


designed_mean_RSA_values_array_00 = array(designed_mean_RSA_values_00)
designed_mean_entropy_values_array_00 = array(designed_mean_entropy_values_00)
designed_cor_entropy_RSA_values_array_00 = array(designed_cor_entropy_RSA_values_00)  
designed_mean_KL_values_array_00 = array(designed_mean_KL_values_00)

buried_mean_KL_values_array_00 = array(buried_mean_KL_values_00)
intermediate_mean_KL_values_array_00 = array(intermediate_mean_KL_values_00)
surface_mean_KL_values_array_00 = array(surface_mean_KL_values_00)

buried_mean_entropy_values_array_00 = array(buried_mean_entropy_values_00)
intermediate_mean_entropy_values_array_00 = array(intermediate_mean_entropy_values_00)
surface_mean_entropy_values_array_00 = array(surface_mean_entropy_values_00)

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_0.0_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_00[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_00[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_00[length_counter] 
    designed_mean_KL = designed_mean_KL_values_00[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_0.0_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_entropy\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_entropy = buried_mean_entropy_values_00[length_counter]
    designed_mean_KL = buried_mean_KL_values_00[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy)+ "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

designed_mean_RSA_values_array_03 = array(designed_mean_RSA_values_03)
designed_mean_entropy_values_array_03 = array(designed_mean_entropy_values_03)
designed_cor_entropy_RSA_values_array_03 = array(designed_cor_entropy_RSA_values_03) 
designed_mean_KL_values_array_03 = array(designed_mean_KL_values_03)

buried_mean_KL_values_array_03 = array(buried_mean_KL_values_03)
intermediate_mean_KL_values_array_03 = array(intermediate_mean_KL_values_03)
surface_mean_KL_values_array_03 = array(surface_mean_KL_values_03)

buried_mean_entropy_values_array_03 = array(buried_mean_entropy_values_03)
intermediate_mean_entropy_values_array_03 = array(intermediate_mean_entropy_values_03)
surface_mean_entropy_values_array_03 = array(surface_mean_entropy_values_03)

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_0.3_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_03[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_03[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_03[length_counter] 
    designed_mean_KL = designed_mean_KL_values_03[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()


designed_mean_RSA_values_array_06 = array(designed_mean_RSA_values_06)
designed_mean_entropy_values_array_06 = array(designed_mean_entropy_values_06)
designed_cor_entropy_RSA_values_array_06 = array(designed_cor_entropy_RSA_values_06)
designed_mean_KL_values_array_06 = array(designed_mean_KL_values_06) 

buried_mean_KL_values_array_06 = array(buried_mean_KL_values_06)
intermediate_mean_KL_values_array_06 = array(intermediate_mean_KL_values_06)
surface_mean_KL_values_array_06 = array(surface_mean_KL_values_06)

buried_mean_entropy_values_array_06 = array(buried_mean_entropy_values_06)
intermediate_mean_entropy_values_array_06 = array(intermediate_mean_entropy_values_06)
surface_mean_entropy_values_array_06 = array(surface_mean_entropy_values_06)

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_0.6_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_06[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_06[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_06[length_counter] 
    designed_mean_KL = designed_mean_KL_values_06[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) +"\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_0.6_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_entropy\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_entropy = buried_mean_entropy_values_06[length_counter]
    designed_mean_KL = buried_mean_KL_values_06[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy)+ "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

designed_mean_RSA_values_array_09 = array(designed_mean_RSA_values_09)
designed_mean_entropy_values_array_09 = array(designed_mean_entropy_values_09)
designed_cor_entropy_RSA_values_array_09 = array(designed_cor_entropy_RSA_values_09) 
designed_mean_KL_values_array_09 = array(designed_mean_KL_values_09)

buried_mean_KL_values_array_09 = array(buried_mean_KL_values_09)
intermediate_mean_KL_values_array_09 = array(intermediate_mean_KL_values_09)
surface_mean_KL_values_array_09 = array(surface_mean_KL_values_09)

buried_mean_entropy_values_array_09 = array(buried_mean_entropy_values_09)
intermediate_mean_entropy_values_array_09 = array(intermediate_mean_entropy_values_09)
surface_mean_entropy_values_array_09 = array(surface_mean_entropy_values_09)

length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_0.9_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_09[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_09[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_09[length_counter] 
    designed_mean_KL = designed_mean_KL_values_09[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) +  "\t" + str(designed_mean_KL)+ "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_0.9_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_entropy\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_entropy = buried_mean_entropy_values_09[length_counter]
    designed_mean_KL = buried_mean_KL_values_09[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy)+ "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

designed_mean_RSA_values_array_12 = array(designed_mean_RSA_values_12)
designed_mean_entropy_values_array_12 = array(designed_mean_entropy_values_12)
designed_cor_entropy_RSA_values_array_12 = array(designed_cor_entropy_RSA_values_12) 
designed_mean_KL_values_array_12 = array(designed_mean_KL_values_12)

buried_mean_KL_values_array_12 = array(buried_mean_KL_values_12)
intermediate_mean_KL_values_array_12 = array(intermediate_mean_KL_values_12)
surface_mean_KL_values_array_12 = array(surface_mean_KL_values_12)

buried_mean_entropy_values_array_12 = array(buried_mean_entropy_values_12)
intermediate_mean_entropy_values_array_12 = array(intermediate_mean_entropy_values_12)
surface_mean_entropy_values_array_12 = array(surface_mean_entropy_values_12)

length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_1.2_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_12[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_12[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_12[length_counter] 
    designed_mean_KL = designed_mean_KL_values_12[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_1.2_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_entropy\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_entropy = buried_mean_entropy_values_12[length_counter]
    designed_mean_KL = buried_mean_KL_values_12[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy)+ "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

designed_mean_RSA_values_array_18 = array(designed_mean_RSA_values_18)
designed_mean_entropy_values_array_18 = array(designed_mean_entropy_values_18)
designed_cor_entropy_RSA_values_array_18 = array(designed_cor_entropy_RSA_values_18) 
designed_mean_KL_values_array_18 = array(designed_mean_KL_values_18)

buried_mean_KL_values_array_18 = array(buried_mean_KL_values_18)
intermediate_mean_KL_values_array_18 = array(intermediate_mean_KL_values_18)
surface_mean_KL_values_array_18 = array(surface_mean_KL_values_18)

buried_mean_entropy_values_array_18 = array(buried_mean_entropy_values_18)
intermediate_mean_entropy_values_array_18 = array(intermediate_mean_entropy_values_18)
surface_mean_entropy_values_array_18 = array(surface_mean_entropy_values_18)

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_1.8_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_18[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_18[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_18[length_counter] 
    designed_mean_KL = designed_mean_KL_values_18[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

designed_mean_RSA_values_array_24 = array(designed_mean_RSA_values_24)
designed_mean_entropy_values_array_24 = array(designed_mean_entropy_values_24)
designed_cor_entropy_RSA_values_array_24 = array(designed_cor_entropy_RSA_values_24) 
designed_mean_KL_values_array_24 = array(designed_mean_KL_values_24)

buried_mean_KL_values_array_24 = array(buried_mean_KL_values_24)
intermediate_mean_KL_values_array_24 = array(intermediate_mean_KL_values_24)
surface_mean_KL_values_array_24 = array(surface_mean_KL_values_24)

buried_mean_entropy_values_array_24 = array(buried_mean_entropy_values_24)
intermediate_mean_entropy_values_array_24 = array(intermediate_mean_entropy_values_24)
surface_mean_entropy_values_array_24 = array(surface_mean_entropy_values_24)

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_2.4_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_24[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_24[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_24[length_counter] 
    designed_mean_KL = designed_mean_KL_values_24[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

designed_mean_RSA_values_array_soft = array(designed_mean_RSA_values_soft)
designed_mean_entropy_values_array_soft = array(designed_mean_entropy_values_soft)
designed_cor_entropy_RSA_values_array_soft = array(designed_cor_entropy_RSA_values_soft) 
designed_mean_KL_values_array_soft = array(designed_mean_KL_values_24)

buried_mean_KL_values_array_soft = array(buried_mean_KL_values_soft)
intermediate_mean_KL_values_array_soft = array(intermediate_mean_KL_values_soft)
surface_mean_KL_values_array_soft = array(surface_mean_KL_values_soft)

buried_mean_entropy_values_array_soft = array(buried_mean_entropy_values_soft)
intermediate_mean_entropy_values_array_soft = array(intermediate_mean_entropy_values_soft)
surface_mean_entropy_values_array_soft = array(surface_mean_entropy_values_soft)

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_soft_noah.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_soft[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_soft[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_soft[length_counter] 
    designed_mean_KL = designed_mean_KL_values_soft[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

mean_KL_temp_file = open("graph_mean_KL_all_temp_data_noah.csv", "w")
mean_KL_temp_file.write("PDB\tmean_KL_temp_0.0\tmean_KL_soft\tmean_KL_temp_0.3\tmean_KL_temp_0.6\tmean_KL_temp_0.9\tmean_KL_temp_1.2\tmean_KL_temp_1.8\tmean_KL_temp_2.4\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(designed_mean_KL_values_00[length_counter]) + "\t" +  str(designed_mean_KL_values_soft[length_counter]) + "\t" + str(designed_mean_KL_values_03[length_counter]) + "\t" +  str(designed_mean_KL_values_06[length_counter]) + "\t" + str(designed_mean_KL_values_09[length_counter]) + "\t" + str(designed_mean_KL_values_12[length_counter]) + "\t" +  str(designed_mean_KL_values_18[length_counter]) + "\t" + str(designed_mean_KL_values_24[length_counter]) + "\t" + str(natural_mean_split_KL_values[length_counter]) + "\n"
    mean_KL_temp_file.write(string)
    length_counter = length_counter + 1
mean_KL_temp_file.close() 

mean_KL_temp_file = open("graph_mean_KL_buried_temp_data_noah.csv", "w")
mean_KL_temp_file.write("PDB\tmean_KL_temp_0.0\tmean_KL_soft\tmean_KL_temp_0.3\tmean_KL_temp_0.6\tmean_KL_temp_0.9\tmean_KL_temp_1.2\tmean_KL_temp_1.8\tmean_KL_temp_2.4\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(buried_mean_KL_values_00[length_counter]) + "\t" + str(buried_mean_KL_values_soft[length_counter]) + "\t" + str(buried_mean_KL_values_03[length_counter]) + "\t" +  str(buried_mean_KL_values_06[length_counter]) + "\t" + str(buried_mean_KL_values_09[length_counter]) + "\t" + str(buried_mean_KL_values_12[length_counter]) + "\t" + str(buried_mean_KL_values_18[length_counter]) + "\t" + str(buried_mean_KL_values_24[length_counter]) + "\t" + str(buried_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_temp_file.write(string)
    length_counter = length_counter + 1
mean_KL_temp_file.close()  

mean_KL_temp_file = open("graph_mean_KL_intermediate_temp_data_noah.csv", "w")
mean_KL_temp_file.write("PDB\tmean_KL_temp_0.0\tmean_KL_soft\tmean_KL_temp_0.3\tmean_KL_temp_0.6\tmean_KL_temp_0.9\tmean_KL_temp_1.2\tmean_KL_temp_1.8\tmean_KL_temp_2.4\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(intermediate_mean_KL_values_00[length_counter]) + "\t" + str(intermediate_mean_KL_values_soft[length_counter]) + "\t" + str(intermediate_mean_KL_values_03[length_counter]) + "\t" +  str(intermediate_mean_KL_values_06[length_counter]) + "\t" + str(intermediate_mean_KL_values_09[length_counter]) + "\t" + str(intermediate_mean_KL_values_12[length_counter]) + "\t" + str(intermediate_mean_KL_values_18[length_counter]) + "\t" + str(intermediate_mean_KL_values_24[length_counter]) + "\t" + str(intermediate_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_temp_file.write(string)
    length_counter = length_counter + 1
mean_KL_temp_file.close()  

mean_KL_temp_file = open("graph_mean_KL_surface_temp_data_noah.csv", "w")
mean_KL_temp_file.write("PDB\tmean_KL_temp_0.0\tmean_KL_soft\tmean_KL_temp_0.3\tmean_KL_temp_0.6\tmean_KL_temp_0.9\tmean_KL_temp_1.2\tmean_KL_temp_1.8\tmean_KL_temp_2.4\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(surface_mean_KL_values_00[length_counter]) + "\t" + str(surface_mean_KL_values_soft[length_counter]) + "\t" + str(surface_mean_KL_values_03[length_counter]) + "\t" +  str(surface_mean_KL_values_06[length_counter]) + "\t" + str(surface_mean_KL_values_09[length_counter]) + "\t" + str(surface_mean_KL_values_12[length_counter]) + "\t" + str(surface_mean_KL_values_18[length_counter]) + "\t" + str(surface_mean_KL_values_24[length_counter]) + "\t" + str(surface_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_temp_file.write(string)
    length_counter = length_counter + 1
mean_KL_temp_file.close()  

mean_KL_temp_file = open("graph_mean_entropy_buried_temp_data_noah.csv", "w")
mean_KL_temp_file.write("PDB\tmean_KL_temp_0.0\tmean_KL_soft\tmean_KL_temp_0.03\tmean_KL_temp_0.1\tmean_KL_temp_0.3\tmean_KL_temp_0.6\tmean_KL_temp_0.9\tmean_KL_temp_1.2\mean_KL_temp_1.8\tmean_KL_temp_2.4\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(buried_mean_entropy_values_00[length_counter]) + "\t" + str(buried_mean_entropy_values_soft[length_counter]) + "\t" + str(buried_mean_entropy_values_03[length_counter]) + "\t" +  str(buried_mean_entropy_values_06[length_counter]) + "\t" + str(buried_mean_entropy_values_09[length_counter]) + "\t" + str(buried_mean_entropy_values_12[length_counter]) + "\t"+ str(buried_mean_entropy_values_18[length_counter]) + "\t" + str(buried_mean_entropy_values_24[length_counter])  + "\t" + str(buried_natural_mean_entropy_values[length_counter]) + "\n"
    mean_KL_temp_file.write(string)
    length_counter = length_counter + 1
mean_KL_temp_file.close() 

mean_KL_temp_file = open("graph_mean_entropy_intermediate_temp_data_noah.csv", "w")
mean_KL_temp_file.write("PDB\tmean_KL_temp_0.0\tmean_KL_soft\tmean_KL_temp_0.1\tmean_KL_temp_0.3\tmean_KL_temp_0.6\tmean_KL_temp_0.9\tmean_KL_temp_1.2\tmean_KL_temp_1.8\tmean_KL_temp_2.4\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(intermediate_mean_entropy_values_00[length_counter]) + "\t" + str(intermediate_mean_entropy_values_soft[length_counter]) + "\t" + str(intermediate_mean_entropy_values_03[length_counter]) + "\t" + str(intermediate_mean_entropy_values_06[length_counter]) + "\t" + str(intermediate_mean_entropy_values_09[length_counter]) + "\t" + str(intermediate_mean_entropy_values_12[length_counter]) + "\t" + str(intermediate_mean_entropy_values_18[length_counter]) + "\t" + str(intermediate_mean_entropy_values_24[length_counter]) + "\t" + str(intermediate_natural_mean_entropy_values[length_counter]) + "\n"
    mean_KL_temp_file.write(string)
    length_counter = length_counter + 1
mean_KL_temp_file.close() 

mean_KL_temp_file = open("graph_mean_entropy_surface_temp_data_noah.csv", "w")
mean_KL_temp_file.write("PDB\tmean_KL_temp_0.0\tmean_KL_soft\tmean_KL_temp_0.3\tmean_KL_temp_0.6\tmean_KL_temp_0.9\tmean_KL_temp_1.2\tmean_KL_temp_1.8\tmean_KL_temp_2.4\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(surface_mean_entropy_values_00[length_counter]) + "\t" + str(surface_mean_entropy_values_soft[length_counter]) + "\t" + str(surface_mean_entropy_values_03[length_counter]) + "\t" +  str(surface_mean_entropy_values_06[length_counter]) + "\t" + str(surface_mean_entropy_values_09[length_counter]) + "\t" + str(surface_mean_entropy_values_12[length_counter]) + "\t" + str(surface_mean_entropy_values_18[length_counter]) + "\t" + str(surface_mean_entropy_values_24[length_counter]) + "\t" + str(surface_natural_mean_entropy_values[length_counter]) + "\n"
    mean_KL_temp_file.write(string)
    length_counter = length_counter + 1
mean_KL_temp_file.close()  

stats_file = open("stats_data_noah.csv", "w")
stats_file.write("PDB\tChain\tmean_entropy_temp_0.0\tmean_entropy_soft\tmean_entropy_temp_0.3\tmean_entropy_temp_0.6\tmean_entropy_temp_0.9\tmean_entropy_temp_1.2\tmean_entropy_temp_1.8\tmean_entropy_temp_2.4\tmean_entropy_natural\tcor_entropy_RSA_03\tnatural_cor_entropy_RSA_value\n")
length_counter = 0
while(length_counter<designed_name_length):
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy_values_00[length_counter]) + "\t" + str(designed_mean_entropy_values_soft[length_counter]) + "\t" + str(designed_mean_entropy_values_03[length_counter]) +  "\t" + str(designed_mean_entropy_values_06[length_counter]) +  "\t" + str(designed_mean_entropy_values_09[length_counter]) +  "\t" + str(designed_mean_entropy_values_12[length_counter]) +  "\t" + str(designed_mean_entropy_values_18[length_counter]) +  "\t" + str(designed_mean_entropy_values_24[length_counter]) + "\t"  +str(natural_mean_entropy_values[length_counter]) + "\t" + str(designed_cor_entropy_RSA_values_03[length_counter]) + "\t"+ str(natural_cor_entropy_RSA_values[length_counter]) + "\n"
    stats_file.write(designed_filestring)
    length_counter = length_counter + 1
stats_file.close()



