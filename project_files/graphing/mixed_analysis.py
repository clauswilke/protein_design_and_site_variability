import os, re, string, math
import numpy as np
import analysis_functions as af 
from pylab import *
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

#Import PDB dict
PDBS = ["1B4T_A", "1CI0_A", "1EFV_B", "1G58_B", "1GV3_A", "1HUJ_A", "1HUR_A", "1IBS_A", "1JLW_A", "1KY2_A", "1KZL_A", "1M3U_A", "1MOZ_A", "1OKC_A", "1PV1_A", "1QMV_A", "1R6M_A", "1RII_A", "1V9S_B", "1W7W_B", "1X1O_B", "1XTD_A", "1YPI_A", "1YSB_A", "1ZNN_A", "1ZWK_A", "2A84_A", "2AIU_A", "2BCG_Y", "2BR9_A", "2CFE_A", "2CJM_C", "2CNV_A", "2ESF_A", "2EU8_A", "2FLI_A", "2G0N_B", "2GV5_D"]

def get_mixed_entropy_values(PDB, buried_temp, surface_temp):
    new_entropies = []
    #Make the files
    buried_file  = "align_data_array_" + PDB + "_" + str(buried_temp) +  ".dat"
    surface_file  = "align_data_array_" + PDB + "_" + str(surface_temp) +  ".dat"
    #Get the RSA Values
    RSA = af.make_array(af.get_RSA_Values(buried_file))
    buried_entropies = af.get_native_entropy(buried_file)
    surface_entropies = af.get_native_entropy(surface_file)
    #Get the entropy values
    for i in xrange(len(RSA)):
        if (float(RSA[i]) <=0.25):
            new_entropies.append(buried_entropies[i])
        else:
            new_entropies.append(surface_entropies[i])
    return RSA, new_entropies
                               
                       
            
cor_values1 = []
cor_values2 = []
natural_cor_values = []
rcParams['font.size'] = 20      
rcParams['figure.figsize'] = [11,8]

for PDB in PDBS:
    RSA1, entropy_mix1 = get_mixed_entropy_values(PDB, 0.3, 0.6)
    RSA2, entropy_mix2 = get_mixed_entropy_values(PDB, 0.3, 0.9)
 
    [cor_entropy_RSA_mix1, pvalue1] = pearsonr(RSA1, entropy_mix1)
    cor_entropy_RSA_mix1 = float(cor_entropy_RSA_mix1)
    cor_values1.append(cor_entropy_RSA_mix1)
    
    [cor_entropy_RSA_mix2, pvalue2] = pearsonr(RSA2, entropy_mix2)
    cor_entropy_RSA_mix2 = float(cor_entropy_RSA_mix2)
    cor_values2.append(cor_entropy_RSA_mix2)
    
    natural_file  = "align_natural_data_array_" + PDB + ".dat"
    natural_RSA = af.make_array(af.get_RSA_Values(natural_file))
    natural_entropy = af.get_native_entropy(natural_file)
    [natural_cor_entropy_RSA, pvalue3] = pearsonr(natural_RSA, natural_entropy)
    natural_cor_entropy_RSA = float(natural_cor_entropy_RSA)
    natural_cor_values.append(natural_cor_entropy_RSA)
    

#print cor_values1, cor_values2, natural_cor_values
#This makes the combo boxplot and line plot for the correlation between RSA and Entropy

fig = plt.figure(1, dpi = 400) #, figsize = (14,6))

#rcParams['lines.linewidth'] = 2
ax = axes([0.09, 0.115, 0.85, 0.85])
#text(-0.37, 0.6, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
b1 = boxplot([cor_values1, cor_values2, natural_cor_values], sym = 'ko')
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
xlabel("Temperature")
ylabel("RSA - Entropy Correlation")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
plt.xticks([1, 2, 3], ["0.3, 0.6", "0.3, 0.9", "NS"])
plt.savefig("Duncan_Mixed_Temp_Correlation_Plot.pdf", format = None)
#plt.show()
