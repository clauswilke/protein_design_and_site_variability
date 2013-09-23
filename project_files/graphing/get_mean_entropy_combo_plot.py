import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions
import matplotlib.cm as cm

all_temp_entropy_values = []
all_temp_entropy_values_noah = []
count = 0
natural_data = []
pdb_names_noah = []
chain_names_noah = []


natural_mean_entropy_values = []
natural_mean_entropy_values_noah = []
designed_mean_entropy_values_00 = []
designed_mean_entropy_values_01 = []
designed_mean_entropy_values_03 = []
designed_mean_entropy_values_06 = []
designed_mean_entropy_values_09 = []
designed_mean_entropy_values_12 = []
designed_mean_entropy_values_003 = []


#This saves the figure 
def save_figure(figure_title, PDB):
    save_fig_title = figure_title + PDB +  ".pdf"
    savefig(save_fig_title, format = None)

protein_file = open("graph_mean_data_natural_noah.csv", "r")
natural_protein_data = protein_file.readlines()
protein_file.close()

header = natural_protein_data.pop(0)
for line in natural_protein_data:
    data = re.split("\t", line)
    natural_data.append(data)

for data in natural_data:
    #pdb_names_noah.append(data[0])
    #chain_names_noah.append(data[1])
    #natural_mean_RSA_values_noah.append(data[2])
    natural_mean_entropy_values_noah.append(data[3])
    #natural_cor_entropy_RSA_values_noah.append(data[4])
    #natural_mean_split_KL_values_noah.append(data[5])

#natural_mean_entropy_values_array = analysis_functions.make_array(natural_mean_entropy_values)
#print natural_mean_entropy_values_noah

protein_file_name = "graph_mean_data_0.0_noah.csv"
[designed_mean_RSA_values_00, designed_mean_entropy_values_00, designed_cor_entropy_RSA_values_00, designed_mean_KL_values_00] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_00 = array(designed_mean_entropy_values_00)


protein_file_name = "graph_mean_data_0.3_noah.csv"
[designed_mean_RSA_values_03, designed_mean_entropy_values_03, designed_cor_entropy_RSA_values_03, designed_mean_KL_values_03] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_03 = array(designed_mean_entropy_values_03)


protein_file_name = "graph_mean_data_0.6_noah.csv"
[designed_mean_RSA_values_06, designed_mean_entropy_values_06, designed_cor_entropy_RSA_values_06, designed_mean_KL_values_06] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_06 = array(designed_mean_entropy_values_06)
 

protein_file_name = "graph_mean_data_0.9_noah.csv"
[designed_mean_RSA_values_09, designed_mean_entropy_values_09, designed_cor_entropy_RSA_values_09, designed_mean_KL_values_09] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_09 = array(designed_mean_entropy_values_09)

protein_file_name = "graph_mean_data_1.2_noah.csv"
[designed_mean_RSA_values_12, designed_mean_entropy_values_12, designed_cor_entropy_RSA_values_12, designed_mean_KL_values_12] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_12 = array(designed_mean_entropy_values_12)

protein_file_name = "graph_mean_data_1.8_noah.csv"
[designed_mean_RSA_values_18, designed_mean_entropy_values_18, designed_cor_entropy_RSA_values_18, designed_mean_KL_values_18] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_18 = array(designed_mean_entropy_values_18)


protein_file_name = "graph_mean_data_2.4_noah.csv"
[designed_mean_RSA_values_24, designed_mean_entropy_values_24, designed_cor_entropy_RSA_values_24, designed_mean_KL_values_24] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_24 = array(designed_mean_entropy_values_24)


protein_file_name = "graph_mean_data_soft_noah.csv"
[designed_mean_RSA_values_soft, designed_mean_entropy_values_soft, designed_cor_entropy_RSA_values_soft, designed_mean_KL_values_soft] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_soft = array(designed_mean_entropy_values_soft)

all_temp_entropy_values_noah.append(designed_mean_entropy_values_00)
all_temp_entropy_values_noah.append(designed_mean_entropy_values_soft)
all_temp_entropy_values_noah.append(designed_mean_entropy_values_03)
all_temp_entropy_values_noah.append(designed_mean_entropy_values_06)
all_temp_entropy_values_noah.append(designed_mean_entropy_values_09)
all_temp_entropy_values_noah.append(designed_mean_entropy_values_12)
all_temp_entropy_values_noah.append(designed_mean_entropy_values_18)
all_temp_entropy_values_noah.append(designed_mean_entropy_values_24)
all_temp_entropy_values_noah.append(natural_mean_entropy_values_noah)

#print all_temp_entropy_values_noah
#print natural_mean_entropy_values_noah

natural_data = []
pdb_names = []
chain_names = []

protein_file = open("graph_mean_data_natural.csv", "r")
natural_protein_data = protein_file.readlines()
protein_file.close()
#print natural_protein_data
header = natural_protein_data.pop(0)
for line in natural_protein_data:
    data = re.split("\t", line)
    natural_data.append(data)

for data in natural_data:
    pdb_names.append(data[0])
    chain_names.append(data[1])
    natural_mean_entropy_values.append(data[3])

#natural_mean_RSA_values_array = analysis_functions.make_array(natural_mean_RSA_values)
natural_mean_entropy_values_array = analysis_functions.make_array(natural_mean_entropy_values)
#natural_cor_entropy_RSA_values_array = analysis_functions.make_array(natural_cor_entropy_RSA_values)
#natural_mean_split_KL_values_array = analysis_functions.make_array(natural_mean_split_KL_values)

protein_file_name = "graph_mean_data_0.0.csv"
[designed_mean_RSA_values_00, designed_mean_entropy_values_00, designed_cor_entropy_RSA_values_00, designed_mean_KL_values_00] = analysis_functions.get_mean_designed_data(protein_file_name)

#designed_mean_RSA_values_array_00 = array(designed_mean_RSA_values_00)
designed_mean_entropy_values_array_00 = array(designed_mean_entropy_values_00)
#designed_cor_entropy_RSA_values_array_00 = array(designed_cor_entropy_RSA_values_00)  
#designed_mean_KL_values_array_00 = array(designed_mean_KL_values_00)

protein_file_name = "graph_mean_data_0.03.csv"
[designed_mean_RSA_values_003, designed_mean_entropy_values_003, designed_cor_entropy_RSA_values_003, designed_mean_KL_values_003] = analysis_functions.get_mean_designed_data(protein_file_name)
#print natural_mean_entropy_values
#designed_mean_RSA_values_array_003 = array(designed_mean_RSA_values_003)
designed_mean_entropy_values_array_003 = array(designed_mean_entropy_values_003)
#designed_cor_entropy_RSA_values_array_003 = array(designed_cor_entropy_RSA_values_003) 
#designed_mean_KL_values_array_003 = array(designed_mean_KL_values_003)

protein_file_name = "graph_mean_data_0.1.csv"
[designed_mean_RSA_values_01, designed_mean_entropy_values_01, designed_cor_entropy_RSA_values_01, designed_mean_KL_values_01] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_01 = array(designed_mean_entropy_values_01)

protein_file_name = "graph_mean_data_0.3.csv"
[designed_mean_RSA_values_03, designed_mean_entropy_values_03, designed_cor_entropy_RSA_values_03, designed_mean_KL_values_03] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_03 = array(designed_mean_entropy_values_03)


protein_file_name = "graph_mean_data_0.6.csv"
[designed_mean_RSA_values_06, designed_mean_entropy_values_06, designed_cor_entropy_RSA_values_06, designed_mean_KL_values_06] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_06 = array(designed_mean_entropy_values_06)

protein_file_name = "graph_mean_data_0.9.csv"
[designed_mean_RSA_values_09, designed_mean_entropy_values_09, designed_cor_entropy_RSA_values_09, designed_mean_KL_values_09] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_09 = array(designed_mean_entropy_values_09)


protein_file_name = "graph_mean_data_1.2.csv"
[designed_mean_RSA_values_12, designed_mean_entropy_values_12, designed_cor_entropy_RSA_values_12, designed_mean_KL_values_12] = analysis_functions.get_mean_designed_data(protein_file_name)
designed_mean_entropy_values_array_12 = array(designed_mean_entropy_values_12)

all_temp_entropy_values.append(designed_mean_entropy_values_00)
all_temp_entropy_values.append(designed_mean_entropy_values_003)
all_temp_entropy_values.append(designed_mean_entropy_values_01)
all_temp_entropy_values.append(designed_mean_entropy_values_03)
all_temp_entropy_values.append(designed_mean_entropy_values_06)
all_temp_entropy_values.append(designed_mean_entropy_values_09)
all_temp_entropy_values.append(designed_mean_entropy_values_12)
all_temp_entropy_values.append(natural_mean_entropy_values)

#print designed_mean_entropy_values_00
#print natural_mean_entropy_values

all_temp_entropy_values_array = []
all_temp_entropy_values_noah_array = []


for element in all_temp_entropy_values:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_temp_entropy_values_array.append(new_array)

for element in all_temp_entropy_values_noah:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_temp_entropy_values_noah_array.append(new_array)

rcParams['figure.figsize'] = [16,6]
#Make boxplot for Temp vs Mean Entropy 
fig = plt.figure(count, dpi = 500)#, figsize = (11,8))
rcParams['font.size'] = 20     
#rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2
#ax = fig.add_subplot(111)
ax = axes([0.063, 0.115, 0.43, 0.85])
b1 = ax.boxplot(all_temp_entropy_values_array, sym = 'ko')
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
xlabel("Temperature")
ylabel("Mean Entropy")
text(-0.5, 2.75, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
#plt.ylim(0.6, 2.4)
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4], ["0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4"])
plt.ylim(0.15, 2.75)
plt.yticks([0.25, 0.75, 1.25,1.75, 2.25, 2.75], ["0.25", "0.75", "1.25", "1.75", "2.25", "2.75"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])

#Make boxplot for Temp vs Mean Entropy for Noah Dataset    
#fig = plt.figure(count, dpi = 500, figsize = (11,8))
rcParams['font.size'] = 20     

rcParams['lines.linewidth'] = 2
#ax = fig.add_subplot(111)

ax2 = axes([0.563, 0.115, 0.43, 0.85])
b2 = ax2.boxplot(all_temp_entropy_values_noah_array, sym = 'ko')
setp(b2['whiskers'], color = 'black', linestyle = '-')
setp(b2['boxes'], color =  'black')
setp(b2['caps'], color = 'black')
setp(b2['medians'], color = 'black')
setp(b2['fliers'], color  = 'black')
xlabel("Temperature")
#ylabel("Mean Entropy")
text(-0.6, 2.75, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
#plt.ylim(0.6, 2.6)
#plt.yticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6], ["0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4", "2.6"])
#plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])
plt.ylim(0.15, 2.75)
#plt.xlim(0.8, 9.2)
plt.yticks([0.25, 0.75, 1.25,1.75, 2.25, 2.75], ["0.25", "0.75", "1.25", "1.75", "2.25", "2.75"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])

save_fig_title = "Mean_Entropy_vs_Temp_Combo_Boxplot" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1

#plt.show()