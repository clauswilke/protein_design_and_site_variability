import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
from scipy.stats import *
import matplotlib.pyplot as plt
import analysis_functions

natural_data = []
pdb_names = []
chain_names = []

all_temp_cor_entropy_RSA_values = []
all_temp_entropy_values = []
natural_mean_split_KL_values = []

temps = [0.0, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2]
temp_array = array(temps)

noah_modified_temps = [-0.2, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2, 1.4]
noah_modified_temp_array = array(noah_modified_temps)
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

noah_buried_all_temp_mean_entropy_data_array = [] 
noah_intermediate_all_temp_mean_entropy_data_array = [] 
noah_surface_all_temp_mean_entropy_data_array = [] 

noah_buried_medians = []
noah_intermediate_medians = []
noah_surface_medians = []

for line in buried_mean_entropy_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    buried_all_temp_data.append(data_array)
    #mean_entropy_temp_values_array = analysis_functions.make_array(data)
noah_buried_all_temp_mean_entropy_data_array = array(buried_all_temp_data)
print len(noah_buried_all_temp_mean_entropy_data_array)
noah_buried_medians = median(noah_buried_all_temp_mean_entropy_data_array, axis = 0)
print len(noah_buried_medians)

for line in intermediate_mean_entropy_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    intermediate_all_temp_data.append(data_array)
    #mean_entropy_temp_values_array = analysis_functions.make_array(data)
noah_intermediate_all_temp_mean_entropy_data_array = array(intermediate_all_temp_data)
noah_intermediate_medians = median(noah_intermediate_all_temp_mean_entropy_data_array, axis = 0)

for line in surface_mean_entropy_temp_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    surface_all_temp_data.append(data_array)
    #mean_entropy_temp_values_array = analysis_functions.make_array(data)
noah_surface_all_temp_mean_entropy_data_array = array(surface_all_temp_data)
noah_surface_medians = median(noah_surface_all_temp_mean_entropy_data_array, axis = 0)

for element in all_temp_entropy_values:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_temp_entropy_values_array.append(new_array)
    
mean_entropy_temp_file = open("graph_mean_entropy_buried_temp_data.csv", "r")
buried_mean_entropy_temp_data = mean_entropy_temp_file.readlines()
mean_entropy_temp_file.close()
buried_header = buried_mean_entropy_temp_data.pop(0)

mean_entropy_temp_file = open("graph_mean_entropy_intermediate_temp_data.csv", "r")
intermediate_mean_entropy_temp_data = mean_entropy_temp_file.readlines()
mean_entropy_temp_file.close()
intermediate_header = intermediate_mean_entropy_temp_data.pop(0)

mean_entropy_temp_file = open("graph_mean_entropy_surface_temp_data.csv", "r")
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

line_values = empty(9)
count = 0 #Make boxplot for Temp vs Mean Entropy for Buried Residues     
fig = plt.figure(count, dpi = 500, figsize = (16,6))
rcParams['font.size'] = 20     
rcParams['lines.markersize'] = 8
#rcParams['figure.figsize'] = [16,5]
rcParams['lines.linewidth'] = 2
#ax = axes([0.08, 0.12, 0.90, 0.83])
ax = axes([0.063, 0.115, 0.43, 0.85])

p1 = ax.plot([1,2,3,4,5,6,7], buried_medians[0:7], color = 'black', marker = "o", linestyle = '-')
line_values.fill(buried_medians[7])
p2 = ax.plot([0,1,2,3,4,5,6,7,8], line_values, color = 'black', linestyle = '--')
p3 = ax.plot([8],buried_medians[7] , color = 'black', marker = 'o')

p3 = ax.plot([1,2,3,4,5,6,7], intermediate_medians[0:7], color = 'blue', marker = "o", linestyle = '-')
line_values.fill(intermediate_medians[7])
p4 = ax.plot([0,1,2,3,4,5,6,7,8], line_values, color = 'blue', linestyle = '--')
p3 = ax.plot([8],intermediate_medians[7] , color = 'blue', marker = 'o')

p5 = ax.plot([1,2,3,4,5,6,7], surface_medians[0:7], color = 'red', marker = "o", linestyle = '-')
line_values.fill(surface_medians[7])
p6 = ax.plot([0,1,2,3,4,5,6,7,8], line_values, color = 'red', linestyle = '--')
p3 = ax.plot([8],surface_medians[7] , color = 'red', marker = 'o')

l1 = ax.legend([p1, p2, p3],["Buried", "Partially Buried", "Exposed"], "lower right", frameon = False, numpoints = 1, prop = {'size': 15}) #Only one dot in the legend, no frame so "False", two columns so ncol = 2, and the font size is 15.

xlabel("Temperature")
ylabel("Median of Mean Entropy")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.ylim(0.25, 2.5)
plt.xlim(0.8,8.2)
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4], ["0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4"])
plt.yticks([0.25,0.5, 0.75, 1.0, 1.25, 1.5,1.75, 2.0, 2.25, 2.5], ["0.25","0.5", "0.75", "1.0", "1.25", "1.5", "1.75", "2.0", "2.25", "2.5"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
text(0.0, 2.5, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
#save_fig_title = "Mean_Entropy_Position_Lineplot" + ".pdf"
#savefig(save_fig_title, format = None)
count = count + 1
#plt.show()


line_values = empty(10)
#print len(noah_buried_medians)
count = 0 #Make boxplot for Temp vs Mean Entropy for Buried Residues     
#fig = plt.figure(count, dpi = 500, figsize = (11,8))
#rcParams['font.size'] = 20     
#rcParams['figure.figsize'] = [11,8]
#rcParams['lines.linewidth'] = 2
#ax2 = axes([0.10, 0.12, 0.90, 0.83])
ax2 = axes([0.563, 0.115, 0.43, 0.85])
#b4 = ax.boxplot(all_temp_mean_entropy_data_array, sym = 'k-)

p1 = ax2.plot([1,2,3,4,5,6,7,8], noah_buried_medians[0:8], color = 'black', marker = "o", linestyle = '-')
line_values.fill(noah_buried_medians[8])
p2 = ax2.plot([0,1,2,3,4,5,6,7,8,9], line_values, color = 'black', linestyle = '--')
p3 = ax2.plot([9],noah_buried_medians[8] , color = 'black', marker = 'o')


p4 = ax2.plot([1,2,3,4,5,6,7,8], noah_intermediate_medians[0:8], color = 'blue', marker = "o", linestyle = '-')
line_values.fill(noah_intermediate_medians[8])
p5 = ax2.plot([0,1,2,3,4,5,6,7,8,9], line_values, color = 'blue', linestyle = '--')
p6 = ax2.plot([9],noah_intermediate_medians[8] , color = 'blue', marker = 'o')


p7 = ax2.plot([1,2,3,4,5,6,7,8], noah_surface_medians[0:8], color = 'red', marker = "o", linestyle = '-')
line_values.fill(noah_surface_medians[8])
p8 = ax2.plot([0,1,2,3,4,5,6,7,8,9], line_values, color = 'red', linestyle = '--')
p9 = ax2.plot([9],noah_surface_medians[8] , color = 'red', marker = 'o')

l1 = ax.legend([p1, p2, p3],["Buried", "Partially Buried", "Exposed"], "lower right", frameon = False, numpoints = 1, prop = {'size': 15}) #Only one dot in the legend, no frame so "False", two columns so ncol = 2, and the font size is 15.

xlabel("Temperature")
#ylabel("Median of Mean Entropy")
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
plt.ylim(0.25, 2.5)
plt.xlim(0.8,9.2)
plt.yticks([0.25, 0.5, 0.75, 1.0, 1.25, 1.5,1.75, 2.0, 2.25, 2.5], ["0.25","0.5", "0.75", "1.0", "1.25", "1.5", "1.75", "2.0", "2.25", "2.5"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])
text(-0.10, 2.5, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
save_fig_title = "Mean_Entropy_Position_Lineplot_Combo" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#plt.show()


