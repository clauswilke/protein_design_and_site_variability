import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions

#Date Last Updated: June 24, 2013
#Description: This is one of the scripts that graphs the data for the paper.

#matplotlib.rcParams['backend'] = "Qt4Agg"
'''def get_plot_color(color_increment, natural_cor_entropy_RSA_value):
    #Start Blue and Go to Red tuple = RGB
    #Start Blue
    #color_start = (1,0,0) 
    green = 0
    #switch_to_red = int(num_total_plots/2)
    if(natural_cor_entropy_RSA_value < 0.2): 
        red = .95 - (0.001 * natural_cor_entropy_RSA_value)
        blue = 0 
        color = (red, green, blue)
    elif(natural_cor_entropy_RSA_value>=0.2 and natural_cor_entropy_RSA_value<0.4):
        red = 1 - (1.2*natural_cor_entropy_RSA_value)
        blue = 1 - (1.2*natural_cor_entropy_RSA_value)
        color = (red,green,blue)
    elif(natural_cor_entropy_RSA_value>=0.4 and natural_cor_entropy_RSA_value<=1):
        red = .4 - (0.65 * natural_cor_entropy_RSA_value)
        blue = 1 - (0.6*natural_cor_entropy_RSA_value)
        color = (red, green, blue)
    else:
        print "Correlation over 1!!!!!"
    return color
'''

def get_plot_color(natural_cor_entropy_RSA_value):
    if(natural_cor_entropy_RSA_value < 0.1): 
        color = cm.hot_r(natural_cor_entropy_RSA_value + 0.2)
    else: 
        color = cm.hot_r(natural_cor_entropy_RSA_value)
    return color


def get_plot_format(xaxis_label, yaxis_label, ax):
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlim(-0.2, 1.4)
    plt.ylim(0.5, 3)
    plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    plt.yticks([-0.2, 0.0, 0.2, 0.4, 0.6])
    #plt.text(0.9, 0.9, "Text",horizontalalignment = 'center', verticalalignment = 'center', transform = ax.transAxes) #Independent of size since it is relative to the axes

def get_plot_format_cor_plot(xaxis_label, yaxis_label, ax):
 
    xlabel(xaxis_label)
    #ylabel(yaxis_label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlim(-0.3, 1.65)
    plt.ylim(-0.3, 0.6)
    #plt.xticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], ["FB", "0.0","0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "NS"])
    plt.xticks([-0.2, 0.05, 0.3, 0.55, 0.8, 1.05, 1.3, 1.55], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
    plt.yticks([-0.2, 0.0, 0.2, 0.4, 0.6])


def save_figure(figure_title, PDB):
    save_fig_title = figure_title + PDB +  ".pdf"
    savefig(save_fig_title, format = None)

count = 1
natural_data = []

pdb_names = []
chain_names = []
natural_mean_RSA_values = []
natural_mean_entropy_values = []
natural_cor_entropy_RSA_values = []
natural_mean_RSA_values = []

designed_mean_RSA_values_00 = []
designed_mean_RSA_values_01 = []
designed_mean_RSA_values_03 = []
designed_mean_RSA_values_06 = []
designed_mean_RSA_values_09 = []
designed_mean_RSA_values_12 = []
designed_mean_RSA_values_003 = []

designed_mean_entropy_values_00 = []
designed_mean_entropy_values_01 = []
designed_mean_entropy_values_03 = []
designed_mean_entropy_values_06 = []
designed_mean_entropy_values_09 = []
designed_mean_entropy_values_12 = []
designed_mean_entropy_values_003 = []

designed_cor_entropy_RSA_values_00 = []
designed_cor_entropy_RSA_values_01 = []
designed_cor_entropy_RSA_values_03 = []
designed_cor_entropy_RSA_values_06 = []
designed_cor_entropy_RSA_values_09 = []
designed_cor_entropy_RSA_values_12 = []
designed_cor_entropy_RSA_values_003 = []

designed_mean_KL_values_00 = []
designed_mean_KL_values_003 = []
designed_mean_KL_values_01 = []
designed_mean_KL_values_03 = []
designed_mean_KL_values_06 = []
designed_mean_KL_values_09 = []
designed_mean_KL_values_12 = []


all_temp_cor_entropy_RSA_values = []
all_temp_entropy_values = []
natural_mean_split_KL_values = []

temps = [0.0, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2]
temp_array = array(temps)

#modified_temps = [-0.2, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2, 1.4]
modified_temps = [-0.2, 0.05, 0.3, 0.55, 0.8, 1.05, 1.3, 1.55]
modified_temp_array = array(modified_temps)

#modified_temps = [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
#modified_temp_array = array(modified_temps)

protein_file = open("graph_mean_data_natural.csv", "r")
natural_protein_data = protein_file.readlines()
protein_file.close()

header = natural_protein_data.pop(0)
for line in natural_protein_data:
    data = re.split("\t", line)
    natural_data.append(data)

for data in natural_data:
    pdb_names.append(data[0])
    chain_names.append(data[1])
    natural_mean_RSA_values.append(data[2])
    natural_mean_entropy_values.append(data[3])
    natural_cor_entropy_RSA_values.append(data[4])
    natural_mean_split_KL_values.append(data[5])

natural_mean_RSA_values_array = analysis_functions.make_array(natural_mean_RSA_values)
natural_mean_entropy_values_array = analysis_functions.make_array(natural_mean_entropy_values)
natural_cor_entropy_RSA_values_array = analysis_functions.make_array(natural_cor_entropy_RSA_values)
natural_mean_split_KL_values_array = analysis_functions.make_array(natural_mean_split_KL_values)

protein_file_name = "graph_mean_data_0.0.csv"
[designed_mean_RSA_values_00, designed_mean_entropy_values_00, designed_cor_entropy_RSA_values_00, designed_mean_KL_values_00] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_00 = array(designed_mean_RSA_values_00)
designed_mean_entropy_values_array_00 = array(designed_mean_entropy_values_00)
designed_cor_entropy_RSA_values_array_00 = array(designed_cor_entropy_RSA_values_00)  
designed_mean_KL_values_array_00 = array(designed_mean_KL_values_00)

protein_file_name = "graph_mean_data_0.03.csv"
[designed_mean_RSA_values_003, designed_mean_entropy_values_003, designed_cor_entropy_RSA_values_003, designed_mean_KL_values_003] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_003 = array(designed_mean_RSA_values_003)
designed_mean_entropy_values_array_003 = array(designed_mean_entropy_values_003)
designed_cor_entropy_RSA_values_array_003 = array(designed_cor_entropy_RSA_values_003) 
designed_mean_KL_values_array_003 = array(designed_mean_KL_values_003)

protein_file_name = "graph_mean_data_0.1.csv"
[designed_mean_RSA_values_01, designed_mean_entropy_values_01, designed_cor_entropy_RSA_values_01, designed_mean_KL_values_01] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_01 = array(designed_mean_RSA_values_01)
designed_mean_entropy_values_array_01 = array(designed_mean_entropy_values_01)
designed_cor_entropy_RSA_values_array_01 = array(designed_cor_entropy_RSA_values_01)  
designed_mean_KL_values_array_01 = array(designed_mean_KL_values_01)

protein_file_name = "graph_mean_data_0.3.csv"
[designed_mean_RSA_values_03, designed_mean_entropy_values_03, designed_cor_entropy_RSA_values_03, designed_mean_KL_values_03] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_03 = array(designed_mean_RSA_values_03)
designed_mean_entropy_values_array_03 = array(designed_mean_entropy_values_03)
designed_cor_entropy_RSA_values_array_03 = array(designed_cor_entropy_RSA_values_03) 
designed_mean_KL_values_array_03 = array(designed_mean_KL_values_03)

protein_file_name = "graph_mean_data_0.6.csv"
[designed_mean_RSA_values_06, designed_mean_entropy_values_06, designed_cor_entropy_RSA_values_06, designed_mean_KL_values_06] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_06 = array(designed_mean_RSA_values_06)
designed_mean_entropy_values_array_06 = array(designed_mean_entropy_values_06)
designed_cor_entropy_RSA_values_array_06 = array(designed_cor_entropy_RSA_values_06)
designed_mean_KL_values_array_06 = array(designed_mean_KL_values_06) 

protein_file_name = "graph_mean_data_0.9.csv"
[designed_mean_RSA_values_09, designed_mean_entropy_values_09, designed_cor_entropy_RSA_values_09, designed_mean_KL_values_09] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_09 = array(designed_mean_RSA_values_09)
designed_mean_entropy_values_array_09 = array(designed_mean_entropy_values_09)
designed_cor_entropy_RSA_values_array_09 = array(designed_cor_entropy_RSA_values_09) 
designed_mean_KL_values_array_09 = array(designed_mean_KL_values_09)

protein_file_name = "graph_mean_data_1.2.csv"
[designed_mean_RSA_values_12, designed_mean_entropy_values_12, designed_cor_entropy_RSA_values_12, designed_mean_KL_values_12] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_12 = array(designed_mean_RSA_values_12)
designed_mean_entropy_values_array_12 = array(designed_mean_entropy_values_12)
designed_cor_entropy_RSA_values_array_12 = array(designed_cor_entropy_RSA_values_12) 
designed_mean_KL_values_array_12 = array(designed_mean_KL_values_12)

all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_00)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_003)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_01)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_03)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_06)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_09)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_12)
all_temp_cor_entropy_RSA_values.append(natural_cor_entropy_RSA_values)

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

all_temp_cor_entropy_RSA_values_length = len(all_temp_cor_entropy_RSA_values[0])
all_temp_cor_entropy_RSA_values_array = []
all_temp_cor_entropy_RSA_values_array = []
all_temp_entropy_values_array = []

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


mean_KL_temp_file = open("graph_mean_KL_all_temp_data.csv", "r")
mean_KL_temp_data = mean_KL_temp_file.readlines()
mean_KL_temp_file.close()
header = mean_KL_temp_data.pop(0)

mean_KL_temp_ordered_file = open("graph_mean_KL_all_temp_data_ordered.csv", "r")
mean_KL_temp_ordered_data = mean_KL_temp_ordered_file.readlines()
mean_KL_temp_ordered_file.close()
ordered_header = mean_KL_temp_ordered_data.pop(0)

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

#This creates the combined boxplot figure for the Mean KL
rcParams['font.size'] = 20     
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2

fig2 = plt.figure(count, dpi = 400)
ax1 = axes([0.07, 0.118, 0.43, 0.849])
text(-0.375, 6.0, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
grid()
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
text(8.75, 6.0, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], ["0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
#count = count + 1

#Make boxplot for Temp vs Mean KL with ordered boxplot  
grid()
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
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], ["0.0", "1.0", "2.0", "3.0","4.0", "5.0", "6.0"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
save_fig_title = "Mean_KL_vs_Temp_Boxplot" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
 
fig = plt.figure(count, dpi = 500)
rcParams['font.size'] = 20     
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2
ax = fig.add_subplot(111)
b3 = ax.boxplot(all_temp_cor_entropy_RSA_values_array, sym = 'ko')
setp(b3['whiskers'], color = 'black', linestyle = '-') #Formats the box boxplots... All Black Everything
setp(b3['boxes'], color =  'black')
setp(b3['caps'], color = 'black')
setp(b3['medians'], color = 'black')
setp(b3['fliers'], color  = 'black')
xlabel("Temperature")
ylabel("RSA - Entropy Correlation")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.get_xaxis().tick_bottom() #Turns off the bottom and left ticks so they are invisible
ax.get_yaxis().tick_left()
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
save_fig_title = "Cor_Mean_Entropy_RSA_vs_Temp_Boxplot" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#show()

#Make boxplot for Temp vs Mean Entropy 
fig = plt.figure(count, dpi = 500, figsize = (11,8))
rcParams['font.size'] = 20     
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2
#ax = fig.add_subplot(111)
ax = axes([0.08,0.12,0.90, 0.83])
b4 = ax.boxplot(all_temp_entropy_values_array, sym = 'ko')
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
save_fig_title = "Mean_Entropy_vs_Temp_Boxplot" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1



all_temp_cor_entropy_RSA_values_length = len(all_temp_cor_entropy_RSA_values[0])
all_temp_cor_entropy_RSA_values_array = []
index = 0

for element in all_temp_cor_entropy_RSA_values:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_temp_cor_entropy_RSA_values_array.append(new_array)

all_temp_cor_entropy_values_transpose = transpose(all_temp_cor_entropy_RSA_values_array)

#This makes the combo boxplot and line plot for the correlation between RSA and Entropy
(m,n) = all_temp_cor_entropy_values_transpose.shape
fig4 = plt.figure(count, dpi = 400) #, figsize = (14,6))
rcParams['font.size'] = 20      
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2
i = 0
ax = axes([0.075, 0.10, 0.43, 0.85])
text(-0.37, 0.6, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
b3 = boxplot(all_temp_cor_entropy_RSA_values_array, sym = 'ko')
setp(b3['whiskers'], color = 'black', linestyle = '-')
setp(b3['boxes'], color =  'black')
setp(b3['caps'], color = 'black')
setp(b3['medians'], color = 'black')
setp(b3['fliers'], color  = 'black')
xlabel("Temperature")
ylabel("RSA - Entropy Correlation")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8], ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
m_first = m - 1 
m_last = m
while (i < m_first): 
    #print i
    PDB = pdb_names[i]
    cor_entropy_RSA_array = all_temp_cor_entropy_values_transpose[i,:]
    natural_cor_entropy_RSA = cor_entropy_RSA_array[7]
    #natural_cor_entropy_RSA[0] = natural_cor_entropy_RSA[0] +  0.05
    #color_tuple = get_plot_color(i, natural_cor_entropy_RSA)
    ax2 = axes([0.565, 0.10, 0.43, 0.85])   
    get_plot_format_cor_plot("Temperature", "RSA - Entropy Correlation", ax2)
    #setp(ax2.get_yticklabels(), visible = False)
    #p1 = plot(modified_temp_array, cor_entropy_RSA_array, color = color_tuple, linestyle = "-", marker = "o")
    p1 = plot(modified_temp_array, cor_entropy_RSA_array, color = get_plot_color(natural_cor_entropy_RSA), linestyle = "-", marker = "o")
    i = i + 1
text(-0.49, 0.6, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
#plt.xticks(modified_temp_array, ["FB", "0.03", "0.1", "0.3", "0.6", "0.9", "1.2", "NS"])
save_fig_title = "Cor_Mean_Entropy_RSA_Combination_Plot" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#plt.show()
