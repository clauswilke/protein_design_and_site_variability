import sys, os, math, string #, re, gzip, urllib, shutil, Bio, subprocess
import analysis_functions as af
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

temps = [0.0, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2]
PDBS = ["1B4T_A", "1CI0_A", "1EFV_B", "1G58_B", "1GV3_A", "1HUJ_A", "1HUR_A", "1IBS_A", "1JLW_A", "1KY2_A", "1KZL_A", "1M3U_A", "1MOZ_A", "1OKC_A", "1PV1_A", "1QMV_A", "1R6M_A", "1RII_A", "1V9S_B", "1W7W_B", "1X1O_B", "1XTD_A", "1YPI_A", "1YSB_A", "1ZNN_A", "1ZWK_A", "2A84_A", "2AIU_A", "2BCG_Y", "2BR9_A", "2CFE_A", "2CJM_C", "2CNV_A", "2ESF_A", "2EU8_A", "2FLI_A", "2G0N_B", "2GV5_D"]

AA = resdict.values()
def get_position_count_data(file):
    buried_counts = []
    surface_counts = []
    buried_sites = []
    surface_sites= []
    input = open(file, "r")
    protein_data = input.readlines()
    input.close()

    RSA = af.get_RSA_Values(file)
    alignment_length = len(RSA)
    site_data = af.get_transformed_data(protein_data)
    #print RSA
    #print site_data
    i = 0
    for site in site_data:
        if(float(RSA[i]) < 0.05):
            buried_sites.append(site)
        else:
            surface_sites.append(site)
        i = i + 1

    buried_array = np.array(buried_sites)
    surface_array = np.array(surface_sites)
    #print buried_array
    #print buried_array
    buried_m, buried_n = buried_array.shape
    buried_total_sum = sum(sum(buried_array))
    
    #print buried_m
    #print buried_n
    surface_m, surface_n = surface_array.shape
    surface_total_sum = sum(sum(surface_array))
    #print surface_total_sum

    j = 0
    while (j < buried_n):
        buried_site_sum = sum(buried_array[:, j])
        buried_counts.append(float(buried_site_sum)/float(buried_total_sum))
        j = j + 1
    
    j = 0
    while (j < surface_n):
        surface_site_sum = sum(surface_array[:, j])
        surface_counts.append(float(surface_site_sum)/float(surface_total_sum))
        j = j + 1
    
    return buried_counts, surface_counts

#pdb_id = "1B4T_A"
#temp = 0.0
#file = "align_data_array_" + pdb_id + "_" + str(temp) + ".dat"  
#[buried_aa, surface_aa] = get_position_count_data(file)
#print buried_aa
#print surface_aa


rcParams['font.size'] = 20
rcParams['figure.figsize'] = [14,6]
#rcParams['lines.markersize'] = 8

N = 20 #The number of bars in each group - should be 20 for the 20 amino acids
index = np.arange(N) # The x locations for the groups 

pp = PdfPages("Designed_Duncan_Frequency_Bar_Graphs.pdf")
count = 0
#plt.show()
for temp in temps:
    num_sites = 0
    
    for pdb_id in PDBS:
        file = "align_data_array_" + pdb_id + "_" + str(temp) + ".dat" 
        print file
        [buried_aa, surface_aa] = get_position_count_data(file)
        buried_array = np.array(buried_aa)
        surface_array = np.array(surface_aa)
        #total_buried_aa = total_buried_aa + buried_aa
        #total_surface_aa = total_surface_aa + surface_aa
        
    
fig = plt.figure(count) 
ax = plt.axes([0.09, 0.12, 0.85, 0.85])
#print index
#print buried_aa
width = 0.40 #The width of the bars
b1 = plt.bar(index, buried_aa, width, color = "red")
b2 = plt.bar(index + width, surface_aa, width, color = "blue") 
ax.text(1.5, 0.285, pdb_id, fontweight = 'bold', ha = 'center', fontsize = 18)#, va = 'center', fontsize = 16)
ax.text(0.6, 0.265,"Temp = " + str(temp), fontweight = 'bold', fontsize = 18)#, ha = 'center', va = 'center', fontsize = 16)
ax.set_xticklabels(AA)
ax.set_xticks(index+width)
ax.set_xlabel("Amino Acid")
ax.set_ylabel("Frequency")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
ax.set_yticklabels(["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"])
plt.legend([b1, b2], ["Buried", "Exposed"], numpoints = 1, frameon = False, loc = 1 , prop = {'size': 15})
pp.savefig(fig)
count = count + 1
#if pdb_id == '1B4T_A' and temp == 0.0:
#    plt.show()
#    quit()

pp.close()  

