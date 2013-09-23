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
#print AA
def get_mean_AA_freqs(file):
    freqs = []
    aa_counts = []
    input = open(file, "r")
    protein_data = input.readlines()
    input.close()
    
    site_data = af.get_transformed_data(protein_data)
    for site in site_data:
        new_site = []
        for element in site:
          new_site.append(float(element))
        freqs.append(new_site)
    freq_array = np.array(freqs)
 
    m,n = freq_array.shape
    total_aa = sum(freq_array)
    m,n = freq_array.shape

    j = 0
    while (j < n):
        aa_freq_sum = sum(freq_array[:, j])
        aa_counts.append(float(aa_freq_sum)/float(total_aa))
        j = j + 1
                       
    return aa_counts

file = "align_data_array_1B4T_A_0.0.dat" 

#data = get_mean_AA_freqs(file)
#print data
rcParams['font.size'] = 20
rcParams['figure.figsize'] = [14,6]
#rcParams['lines.markersize'] = 8

N = 20 #The number of bars in each group - should be 20 for the 20 amino acids
index = np.arange(N) # The x locations for the groups 

count = 0
for temp in temps:
    total_aa = []
    for pdb_id in PDBS:
        file = "align_data_array_" + pdb_id + "_" + str(temp) + ".dat" 
        print file
        aa_freqs = get_mean_AA_freqs(file)
        total_aa.append(aa_freqs)
           
      
    aa_freq_array = np.array(total_aa)
    aa_rows, aa_cols = aa_freq_array.shape
    mean_aa_freqs = []

    for i in xrange(0, aa_cols):
        mean_aa_freqs.append(np.mean(aa_freq_array[:, i]))
 
    fig = plt.figure(count, dpi = 500) 
    ax = plt.axes([0.09, 0.12, 0.85, 0.85])
    width = 0.50 #The width of the bars
    b1 = plt.bar(index, mean_aa_freqs, width, color = "black")
    ax.text(1.5, 0.285, "Temp = " + str(temp), fontweight = 'bold', ha = 'center', fontsize = 18)#, va = 'center', fontsize = 16)
    ax.set_xticklabels(AA)
    ax.set_xticks((index+width/2.))
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Frequency")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
    ax.set_yticklabels(["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"])
    fig_title = "Duncan_Total_Frequency_Bar_Graphs_" + str(temp) + ".pdf"
    plt.savefig(fig_title, format = None)
    count = count + 1

#count = 0
total_aa = []
for pdb_id in PDBS:
    file = "align_natural_data_array_" + pdb_id  + ".dat" 
    print file
    aa_freqs = get_mean_AA_freqs(file)
    total_aa.append(aa_freqs)

aa_freq_array = np.array(total_aa)
aa_rows, aa_cols = aa_freq_array.shape
mean_aa_freqs = []

for i in xrange(0, aa_cols):
    mean_aa_freqs.append(np.mean(aa_freq_array[:, i]))

fig = plt.figure(count, dpi = 500) 
ax = plt.axes([0.09, 0.12, 0.85, 0.85])
width = 0.50 #The width of the bars
b1 = plt.bar(index, mean_aa_freqs, width, color = "black")
ax.text(1.5, 0.285, "Natural", fontweight = 'bold', ha = 'center', fontsize = 18)#, va = 'center', fontsize = 16)
ax.set_xticklabels(AA)
ax.set_xticks((index+width/2.))
ax.set_xlabel("Amino Acid")
ax.set_ylabel("Frequency")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
ax.set_yticklabels(["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"])
fig_title = "Duncan_Total_Frequency_Bar_Graphs_Natural" +  ".pdf"
plt.savefig(fig_title, format = None)
count = count + 1



#plt.show()



