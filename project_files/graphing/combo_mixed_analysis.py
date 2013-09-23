import os, re, string, math
import numpy as np
import analysis_functions as af 
from pylab import *
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

#This is dictionary that maps the natural sequence identifier with its PDB for the Noah Dataset
identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076':'2X1B', \
	'PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37', \
	'PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY', 'PF00355' :'1FQT', \
	'PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS', 'PF00542' :'1CTF', \
	'PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8',  \
	'PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI', 'PF02823' :'1AQT', \
	'PF04002':'2QLC','PF07686' :'2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

def get_noah_mixed_entropy_values(PDB, buried_temp, surface_temp):
    new_entropies = []
    #Make the files
    buried_file  = "align_data_array_" + PDB + "_" +  chain_id  + "_" +  str(buried_temp) +  ".dat"
    surface_file  = "align_data_array_" + PDB + "_" + chain_id  + "_" +  str(surface_temp) +  ".dat"
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
                                                                                    
            
#Import PDB dict
PDBS = ["1B4T_A", "1CI0_A", "1EFV_B", "1G58_B", "1GV3_A", "1HUJ_A", "1HUR_A", "1IBS_A", "1JLW_A", "1KY2_A", "1KZL_A", "1M3U_A", "1MOZ_A", "1OKC_A", "1PV1_A", "1QMV_A", "1R6M_A", "1RII_A", "1V9S_B", "1W7W_B", "1X1O_B", "1XTD_A", "1YPI_A", "1YSB_A", "1ZNN_A", "1ZWK_A", "2A84_A", "2AIU_A", "2BCG_Y", "2BR9_A", "2CFE_A", "2CJM_C", "2CNV_A", "2ESF_A", "2EU8_A", "2FLI_A", "2G0N_B", "2GV5_D"]
                                   
cor_values1 = []
cor_values2 = []
natural_cor_values = []
rcParams['font.size'] = 20      

duncan_file = open("mixed_duncan_stats_file.csv", "w")
duncan_file.write("PDB\tcor_entropy_00_01\tcor_entropy_003_01\tnatural_cor_entropy_value\n")
for PDB in PDBS:
    RSA1, entropy_mix1 = get_mixed_entropy_values(PDB, 0.0, 0.1)
    RSA2, entropy_mix2 = get_mixed_entropy_values(PDB, 0.03, 0.1)
 
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

    duncan_file.write(PDB + "\t" + str(cor_entropy_RSA_mix1) + "\t" + str(cor_entropy_RSA_mix2) + "\t" + str(natural_cor_entropy_RSA) + "\n") 
duncan_file.close()

fig = plt.figure(1, dpi = 400, figsize = (16,6))

#rcParams['lines.linewidth'] = 2
ax = axes([0.066, 0.115, 0.43, 0.85])
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
text(0.2, 0.8, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
plt.yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
plt.xticks([1, 2, 3], ["FB, 0.1", "0.03, 0.1", "NS"])
#plt.savefig("Duncan_Mixed_Temp_Correlation_Plot.pdf", format = None)
#plt.show()

cor_values1 = []
cor_values2 = []
natural_cor_values = []
rcParams['font.size'] = 20   
PDBS = identity_dict.values()
chain_id = 'A'

noah_file = open("mixed_noah_stats_file.csv", "w")
noah_file.write("PDB\tcor_entropy_03_18\tcor_entropy_06_18\tnatural_cor_entropy_value\n")
for PDB in PDBS:
    RSA1, entropy_mix1 = get_noah_mixed_entropy_values(PDB, 0.3, 1.8)
    RSA2, entropy_mix2 = get_noah_mixed_entropy_values(PDB, 0.6, 1.8)
 
    [cor_entropy_RSA_mix1, pvalue1] = pearsonr(RSA1, entropy_mix1)
    cor_entropy_RSA_mix1 = float(cor_entropy_RSA_mix1)
    cor_values1.append(cor_entropy_RSA_mix1)
    
    [cor_entropy_RSA_mix2, pvalue2] = pearsonr(RSA2, entropy_mix2)
    cor_entropy_RSA_mix2 = float(cor_entropy_RSA_mix2)
    cor_values2.append(cor_entropy_RSA_mix2)
    
    natural_file  = "align_natural_data_array_" + PDB + "_" +  chain_id + ".dat"
    natural_RSA = af.make_array(af.get_RSA_Values(natural_file))
    natural_entropy = af.get_native_entropy(natural_file)
    [natural_cor_entropy_RSA, pvalue3] = pearsonr(natural_RSA, natural_entropy)
    natural_cor_entropy_RSA = float(natural_cor_entropy_RSA)
    natural_cor_values.append(natural_cor_entropy_RSA) 

    noah_file.write(PDB + "\t" + str(cor_entropy_RSA_mix1) + "\t" + str(cor_entropy_RSA_mix2) + "\t" + str(natural_cor_entropy_RSA) + "\n") 
noah_file.close() 
    
ax2 = axes([0.563, 0.115, 0.43, 0.85])
#text(-0.37, 0.6, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
b2 = boxplot([cor_values1, cor_values2, natural_cor_values], sym = 'ko')
setp(b2['whiskers'], color = 'black', linestyle = '-')
setp(b2['boxes'], color =  'black')
setp(b2['caps'], color = 'black')
setp(b2['medians'], color = 'black')
setp(b2['fliers'], color  = 'black')
xlabel("Temperature")
#ylabel("RSA - Entropy Correlation")
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
text(0.2, 0.8, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
plt.yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
plt.xticks([1, 2, 3], ["0.3, 1.8", "0.6, 1.8", "NS"])
plt.savefig("Combo_Mixed_Temp_Correlation_Plot.pdf", format = None)
#plt.show()



