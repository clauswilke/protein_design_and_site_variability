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
protein_file = open("graph_mean_data_natural_noah.csv", "r")
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
    
mean_entropy_temp_file = open("graph_mean_entropy_buried_temp_data_noah.csv", "r")
buried_mean_entropy_temp_data = mean_entropy_temp_file.readlines()
mean_entropy_temp_file.close()
buried_header = buried_mean_entropy_temp_data.pop(0)

mean_entropy_temp_file = open("graph_mean_entropy_intermediate_temp_data_noah.csv", "r")
intermediate_mean_entropy_temp_data = mean_entropy_temp_file.readlines()
mean_entropy_temp_file.close()
intermediate_header = intermediate_mean_entropy_temp_data.pop(0)

mean_entropy_temp_file = open("graph_mean_entropy_surface_temp_data_noah.csv", "r")
surface_mean_entropy_temp_data = mean_entropy_temp_file.readlines()
mean_entropy_temp_file.close()
surface_header = surface_mean_entropy_temp_data.pop(0)

all_temp_data = []

buried_all_temp_data = []
intermediate_all_temp_data = []
surface_all_temp_data = []

buried_all_temp_mean_entropy_data_array = [] 
intermediate_all_temp_mean_entropy_data_array = [] 
surface_all_temp_mean_entropy_data_array = [] 

buried_medians = []
intermediate_medians = []
surface_medians = []

for line in buried_mean_entropy_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    buried_all_temp_data.append(data_array)
    #mean_entropy_temp_values_array = analysis_functions.make_array(data)
buried_all_temp_mean_entropy_data_array = array(buried_all_temp_data)
buried_medians = median(buried_all_temp_mean_entropy_data_array, axis = 0)

for line in intermediate_mean_entropy_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    intermediate_all_temp_data.append(data_array)
    #mean_entropy_temp_values_array = analysis_functions.make_array(data)
intermediate_all_temp_mean_entropy_data_array = array(intermediate_all_temp_data)
intermediate_medians = median(intermediate_all_temp_mean_entropy_data_array, axis = 0)

for line in surface_mean_entropy_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    surface_all_temp_data.append(data_array)
    #mean_entropy_temp_values_array = analysis_functions.make_array(data)
surface_all_temp_mean_entropy_data_array = array(surface_all_temp_data)
surface_medians = median(surface_all_temp_mean_entropy_data_array, axis = 0)


#print buried_medians
#m,n = all_temp_mean_entropy_data_array.shape
#print m
#print n

count = 0 #Make boxplot for Temp vs Mean Entropy for Buried Residues     
fig = plt.figure(count, dpi = 500, figsize = (11,8))
rcParams['font.size'] = 20     
rcParams['figure.figsize'] = [11,8]
rcParams['lines.linewidth'] = 2
ax = axes([0.10, 0.12, 0.90, 0.83])

#b4 = ax.boxplot(all_temp_mean_entropy_data_array, sym = 'k-')
p1 = ax.plot([1,2,3,4,5,6,7,8,9], buried_medians, color = 'black', marker = "o", linestyle = '-')
p2 = ax.plot([1,2,3,4,5,6,7,8,9], intermediate_medians, color = 'blue', marker = "o", linestyle = '-')
p3 = ax.plot([1,2,3,4,5,6,7,8,9], surface_medians, color = 'red', marker = "o", linestyle = '-')

l1 = ax.legend([p1, p2, p3],["Buried", "Partially Buried", "Exposed"], "lower right", frameon = False, numpoints = 1, prop = {'size': 15}) #Only one dot in the legend, no frame so "False", two columns so ncol = 2, and the font size is 15.

xlabel("Temperature")
ylabel("Median of Mean Entropy")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.ylim(0.6, 2.4)
plt.xlim(0.8,9.2)
plt.yticks([0.5, 0.75, 1.0, 1.25, 1.5,1.75, 2.0, 2.25, 2.5], ["0.5", "0.75", "1.0", "1.25", "1.5", "1.75", "2.0", "2.25", "2.5"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])
save_fig_title = "Mean_Entropy_Position_Lineplot_Noah" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#plt.show()


