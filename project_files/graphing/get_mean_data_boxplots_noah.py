import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions
import matplotlib.cm as cm

#Date Last Updated: June 24, 2013
#Description: This is one of the scripts that graphs the data for the paper.

#matplotlib.rcParams['backend'] = "Qt4Agg"
def get_plot_color(color_increment, natural_cor_entropy_RSA_value):
    #Start Blue and Go to Red tuple = RGB
    #Start Blue
    #color_start = (1,0,0) 
    green = 0
    #switch_to_red = int(num_total_plots/2)
    if(natural_cor_entropy_RSA_value < 0.2): 
        #print "In the First: "
        red = .95 - (0.001 * natural_cor_entropy_RSA_value)
        blue = 0 
        color = (red, green, blue)
        #print "Next One: "
    elif(natural_cor_entropy_RSA_value>=0.2 and natural_cor_entropy_RSA_value<0.4):
        red = 1 - (1.2*natural_cor_entropy_RSA_value)
        blue = 1 - (1.2*natural_cor_entropy_RSA_value)
        color = (red,green,blue)
    elif(natural_cor_entropy_RSA_value>=0.4 and natural_cor_entropy_RSA_value<=1):
        #print "Third One: "
        red = .5 - (0.65 * natural_cor_entropy_RSA_value)
        blue = 1 - (0.6*natural_cor_entropy_RSA_value)
        color = (red, green, blue)
    else:
        print "Correlation over 1!!!!!"
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
    plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4])
    plt.yticks([-0.2, 0.0, 0.2, 0.4, 0.6])
    #plt.text(0.9, 0.9, "Text",horizontalalignment = 'center', verticalalignment = 'center', transform = ax.transAxes) #Independent of size since it is relative to the axes

def get_plot_format_cor_plot(xaxis_label, yaxis_label, ax):
 
    xlabel(xaxis_label)
    #ylabel(yaxis_label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlim(-0.3, 1.5)
    plt.ylim(-0.5, 0.6)
    plt.xticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])
    plt.yticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8], ["-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8"])

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

designed_mean_RSA_values_00 = []
designed_mean_RSA_values_03 = []
designed_mean_RSA_values_06 = []
designed_mean_RSA_values_09 = []
designed_mean_RSA_values_12 = []
designed_mean_RSA_values_18 = []
designed_mean_RSA_values_24 = []
designed_mean_RSA_values_soft = []

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

designed_mean_KL_values_00 = []
designed_mean_KL_values_03 = []
designed_mean_KL_values_06 = []
designed_mean_KL_values_09 = []
designed_mean_KL_values_12 = []
designed_mean_KL_values_18 = []
designed_mean_KL_values_24 = []
designed_mean_KL_values_soft = []

all_temp_cor_entropy_RSA_values = []
all_temp_entropy_values = []
natural_mean_split_KL_values = []

temps = [0.0, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4]
temp_array = array(temps)

#modified_temps = [-0.2, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4, 2.6, 2.8]
modified_temps = [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
modified_temp_array = array(modified_temps)
protein_file = open("graph_mean_data_natural_noah.csv", "r")
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

protein_file_name = "graph_mean_data_0.0_noah.csv"
[designed_mean_RSA_values_00, designed_mean_entropy_values_00, designed_cor_entropy_RSA_values_00, designed_mean_KL_values_00] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_00 = array(designed_mean_RSA_values_00)
designed_mean_entropy_values_array_00 = array(designed_mean_entropy_values_00)
designed_cor_entropy_RSA_values_array_00 = array(designed_cor_entropy_RSA_values_00)  
designed_mean_KL_values_array_00 = array(designed_mean_KL_values_00)

protein_file_name = "graph_mean_data_0.3_noah.csv"
[designed_mean_RSA_values_03, designed_mean_entropy_values_03, designed_cor_entropy_RSA_values_03, designed_mean_KL_values_03] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_03 = array(designed_mean_RSA_values_03)
designed_mean_entropy_values_array_03 = array(designed_mean_entropy_values_03)
designed_cor_entropy_RSA_values_array_03 = array(designed_cor_entropy_RSA_values_03) 
designed_mean_KL_values_array_03 = array(designed_mean_KL_values_03)

protein_file_name = "graph_mean_data_0.6_noah.csv"
[designed_mean_RSA_values_06, designed_mean_entropy_values_06, designed_cor_entropy_RSA_values_06, designed_mean_KL_values_06] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_06 = array(designed_mean_RSA_values_06)
designed_mean_entropy_values_array_06 = array(designed_mean_entropy_values_06)
designed_cor_entropy_RSA_values_array_06 = array(designed_cor_entropy_RSA_values_06)
designed_mean_KL_values_array_06 = array(designed_mean_KL_values_06) 

protein_file_name = "graph_mean_data_0.9_noah.csv"
[designed_mean_RSA_values_09, designed_mean_entropy_values_09, designed_cor_entropy_RSA_values_09, designed_mean_KL_values_09] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_09 = array(designed_mean_RSA_values_09)
designed_mean_entropy_values_array_09 = array(designed_mean_entropy_values_09)
designed_cor_entropy_RSA_values_array_09 = array(designed_cor_entropy_RSA_values_09) 
designed_mean_KL_values_array_09 = array(designed_mean_KL_values_09)

protein_file_name = "graph_mean_data_1.2_noah.csv"
[designed_mean_RSA_values_12, designed_mean_entropy_values_12, designed_cor_entropy_RSA_values_12, designed_mean_KL_values_12] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_12 = array(designed_mean_RSA_values_12)
designed_mean_entropy_values_array_12 = array(designed_mean_entropy_values_12)
designed_cor_entropy_RSA_values_array_12 = array(designed_cor_entropy_RSA_values_12) 
designed_mean_KL_values_array_12 = array(designed_mean_KL_values_12)

protein_file_name = "graph_mean_data_1.8_noah.csv"
[designed_mean_RSA_values_18, designed_mean_entropy_values_18, designed_cor_entropy_RSA_values_18, designed_mean_KL_values_18] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_18 = array(designed_mean_RSA_values_18)
designed_mean_entropy_values_array_18 = array(designed_mean_entropy_values_18)
designed_cor_entropy_RSA_values_array_18 = array(designed_cor_entropy_RSA_values_18) 
designed_mean_KL_values_array_18 = array(designed_mean_KL_values_18)

protein_file_name = "graph_mean_data_2.4_noah.csv"
[designed_mean_RSA_values_24, designed_mean_entropy_values_24, designed_cor_entropy_RSA_values_24, designed_mean_KL_values_24] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_24 = array(designed_mean_RSA_values_24)
designed_mean_entropy_values_array_24 = array(designed_mean_entropy_values_24)
designed_cor_entropy_RSA_values_array_24 = array(designed_cor_entropy_RSA_values_24)  
designed_mean_KL_values_array_24 = array(designed_mean_KL_values_24)

protein_file_name = "graph_mean_data_soft_noah.csv"
[designed_mean_RSA_values_soft, designed_mean_entropy_values_soft, designed_cor_entropy_RSA_values_soft, designed_mean_KL_values_soft] = analysis_functions.get_mean_designed_data(protein_file_name)

designed_mean_RSA_values_array_soft = array(designed_mean_RSA_values_soft)
designed_mean_entropy_values_array_soft = array(designed_mean_entropy_values_soft)
designed_cor_entropy_RSA_values_array_soft = array(designed_cor_entropy_RSA_values_soft)  
designed_mean_KL_values_array_soft = array(designed_mean_KL_values_soft)


all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_00)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_soft)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_03)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_06)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_09)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_12)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_18)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_24)
all_temp_cor_entropy_RSA_values.append(natural_cor_entropy_RSA_values)

all_temp_entropy_values.append(designed_mean_entropy_values_00)
all_temp_entropy_values.append(designed_mean_entropy_values_soft)
all_temp_entropy_values.append(designed_mean_entropy_values_03)
all_temp_entropy_values.append(designed_mean_entropy_values_06)
all_temp_entropy_values.append(designed_mean_entropy_values_09)
all_temp_entropy_values.append(designed_mean_entropy_values_12)
all_temp_entropy_values.append(designed_mean_entropy_values_18)
all_temp_entropy_values.append(designed_mean_entropy_values_24)
all_temp_entropy_values.append(natural_mean_entropy_values)

#print all_temp_entropy_values
print natural_mean_entropy_values

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


mean_KL_temp_file = open("graph_mean_KL_all_temp_data_noah.csv", "r")
mean_KL_temp_data = mean_KL_temp_file.readlines()
mean_KL_temp_file.close()
header = mean_KL_temp_data.pop(0)

mean_KL_temp_ordered_file = open("graph_mean_KL_all_temp_data_ordered_noah.csv", "r")
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
    #print data
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
text(-0.4, 6.0, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
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
text(9.85, 6.0, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], ["0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])
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
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], ["0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2","1.8", "2.4", "NS"])
save_fig_title = "Mean_KL_vs_Temp_Boxplot_Noah" + ".pdf"
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
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "0.1", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "Soft", "NS"])
save_fig_title = "Cor_Mean_Entropy_RSA_vs_Temp_Boxplot_Noah" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#show()

#Make boxplot for Temp vs Mean KL    
fig = plt.figure(count, dpi = 500, figsize = (11,8))
rcParams['font.size'] = 20     
#rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2
#ax = fig.add_subplot(111)
ax = axes([0.095,0.12,0.90, 0.83])
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
#plt.ylim(0.6, 2.6)
#plt.yticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6], ["0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4", "2.6"])
#plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])
plt.ylim(0.15, 2.75)
#plt.xlim(0.8, 9.2)
plt.yticks([0.25, 0.75, 1.25,1.75, 2.25, 2.75], ["0.25", "0.75", "1.25", "1.75", "2.25", "2.75"])
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])

save_fig_title = "Mean_Entropy_vs_Temp_Boxplot_Noah" + ".pdf"
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
ax = axes([0.075, 0.11, 0.43, 0.85])
text(-0.40, 0.8, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
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
#plt.xlim(-0.3, 1.5)
plt.ylim(-0.5, 0.6)
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9], ["FB", "Soft", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "NS"])
plt.yticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8], ["-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8"])
#print modified_temp_array
#print all_temp_cor_entropy_values_transpose[i,:]

while (i < m): 
    #print i
    PDB = pdb_names[i]
    cor_entropy_RSA_array = all_temp_cor_entropy_values_transpose[i,:]
    natural_cor_entropy_RSA = cor_entropy_RSA_array[8]
    #color_tuple = get_plot_color(i, natural_cor_entropy_RSA)
    ax2 = axes([0.565, 0.11, 0.43, 0.85])   
    get_plot_format_cor_plot("Temperature", "RSA - Entropy Correlation", ax2)
    #setp(ax2.get_yticklabels(), visible = False)
    #p1 = plot(modified_temp_array, cor_entropy_RSA_array, color = color_tuple, linestyle = "-", marker = "o")
    p1 = plot(modified_temp_array, cor_entropy_RSA_array, color = cm.hot_r(natural_cor_entropy_RSA + 0.07), linestyle = "-", marker = "o")
    i = i + 1
text(-0.49, 0.8, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
save_fig_title = "Cor_Mean_Entropy_RSA_Combination_Plot_Noah" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#plt.show()
