import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
from scipy.stats import *
import matplotlib.pyplot as plt
import analysis_functions

#This saves the figure 
def save_figure(figure_title, PDB):
    save_fig_title = figure_title + PDB +  ".pdf"
    savefig(save_fig_title, format = None)

count = 1
natural_data = []

pdb_names = []
chain_names = []

all_temp_cor_entropy_RSA_values = []
all_temp_entropy_values = []
natural_mean_split_KL_values = []


temps = [0.0, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2]
temp_array = array(temps)

modified_temps = [-0.2, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2, 1.4]
modified_temp_array = array(modified_temps)
protein_file = open("graph_mean_data_natural.csv", "r")
natural_protein_data = protein_file.readlines()
protein_file.close()


index = 0

for element in all_temp_cor_entropy_RSA_values:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_temp_cor_entropy_RSA_values_array.append(new_array)

for element in all_temp_entropy_values:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_temp_entropy_values_array.append(new_array)
    
mean_KL_temp_file = open("graph_mean_KL_buried_temp_data.csv", "r")
mean_KL_temp_data = mean_KL_temp_file.readlines()
mean_KL_temp_file.close()
header = mean_KL_temp_data.pop(0)

mean_KL_temp_ordered_file = open("graph_mean_KL_buried_temp_data_ordered.csv", "r")
mean_KL_temp_ordered_data = mean_KL_temp_ordered_file.readlines()
mean_KL_temp_ordered_file.close()
ordered_header = mean_KL_temp_ordered_data.pop(0)

mean_entropy_temp_file = open("graph_mean_entropy_buried_temp_data.csv", "r")
mean_entropy_temp_data = mean_entropy_temp_file.readlines()
mean_entropy_temp_file.close()
header = mean_entropy_temp_data.pop(0)

all_temp_data = []
all_temp_mean_KL_data_array = [] 
for line in mean_KL_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    all_temp_data.append(data_array)
    mean_KL_temp_values_array = analysis_functions.make_array(data)
all_temp_mean_KL_data_array = array(all_temp_data)

all_temp_ordered_data = []
all_temp_ordered_mean_KL_data_array = [] 
for line in mean_KL_temp_ordered_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    all_temp_ordered_data.append(data_array)
    mean_KL_ordered_temp_values_array = analysis_functions.make_array(data)
all_temp_ordered_mean_KL_data_array = array(all_temp_ordered_data)

all_temp_data = []
all_temp_mean_entropy_data_array = [] 
for line in mean_entropy_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    all_temp_data.append(data_array)
    mean_entropy_temp_values_array = analysis_functions.make_array(data)
all_temp_mean_entropy_data_array = array(all_temp_data)

count = 0

#Make boxplot for Temp vs Mean KL for the Unordered and Ordered Plots
rcParams['font.size'] = 20     
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2
fig2 = plt.figure(count, dpi = 400)
ax1 = axes([0.07, 0.118, 0.43, 0.849])
text(-0.375, 3.5, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
#grid()
b1 = boxplot(all_temp_mean_KL_data_array, sym = "ko")
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
xlabel("Temperature")
ylabel("Mean KL Divergence")
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
text(8.75, 3.5, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.ylim(0.0, 3.5)
plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5], ["0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
count = count + 1
#show()

ax2 = axes([0.56, 0.118, 0.43, 0.849])
b2 = boxplot(all_temp_ordered_mean_KL_data_array, sym = "ko")
setp(b2['whiskers'], color = 'black', linestyle = '-')
setp(b2['boxes'], color =  'black')
setp(b2['caps'], color = 'black')
setp(b2['medians'], color = 'black')
setp(b2['fliers'], color  = 'black')
xlabel("Temperature")
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
plt.ylim(0.0, 3.5)
plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5], ["0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
save_fig_title = "Mean_KL_vs_Temp_Buried_Boxplot" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#show()

#Make boxplot for Temp vs Mean Entropy for Buried Residues    
fig = plt.figure(count, dpi = 500, figsize = (11,8))
rcParams['font.size'] = 20     
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2
ax = axes([0.08,0.12,0.90, 0.83])
b4 = ax.boxplot(all_temp_mean_entropy_data_array, sym = 'ko')
setp(b4['whiskers'], color = 'black', linestyle = '-')
setp(b4['boxes'], color =  'black')
setp(b4['caps'], color = 'black')
setp(b4['medians'], color = 'black')
setp(b4['fliers'], color  = 'black')
xlabel("Temperature")
ylabel("Mean Entropy")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.ylim(0.6, 2.4)
plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4], ["0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
save_fig_title = "Mean_Entropy_vs_Temp_Buried_Boxplot" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#plt.show()


