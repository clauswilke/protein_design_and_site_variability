import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
from scipy.stats import *
import matplotlib.pyplot as plt
import analysis_functions

#Last Updated: June 24, 2013
#Description: This a script that graphs the slope figure used in the paper. 
def get_plot_format_intercepts(xaxis_label, yaxis_label, axis):
    xlabel(xaxis_label)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()

    plt.ylim(-1.6, 1.0)
    #plt.xlim(0.1, 1.5)
    plt.xlim(0.1, 1.7)
    #plt.xticks([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], ["0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4"])
    plt.xticks([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6], ["0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6"])
    plt.yticks([-1.5, -1.0, -0.5, 0, 0.5, 1], ["-1.5","-1.0", "-0.5", "0.0", "0.5", "1.0"])

def get_plot_format_entropy(xaxis_label, yaxis_label, ax):
    ax.spines['top'].set_visible(False) #Get rid of the spines on the top and the right. Make them invsible
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom() #Keep the ticks on the left and the right though
    ax.get_yaxis().tick_left()
    xlabel(xaxis_label)
    ylabel(yaxis_label)

    #plt.xlim(0.7,2.4)
    plt.xlim(0.4,2.6)
    plt.ylim(-1.6,1.0)
    #plt.xticks([0.75, 1, 1.25, 1.5, 1.75, 2, 2.25], ["0.75", "1", "1.25", "1.5", "1.75", "2", "2.25"])
    plt.xticks([0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5], ["0.5" ,"0.75", "1", "1.25", "1.5", "1.75", "2", "2.25", "2.5"])
    plt.yticks([-1.5, -1.0, -0.5, 0, 0.5, 1.0], ["-1.5","-1.0", "-0.5", "0.0", "0.5", "1.0"])

def save_figure(figure_title):
    #save_fig_title = figure_title + ".svg"
    #savefig(save_fig_title, format = None)

    save_fig_title = figure_title + ".pdf"
    savefig(save_fig_title, format = None)

count = 1
natural_data = []

pdb_names = []
chain_names = []
natural_mean_RSA_values = []
natural_mean_entropy_values = []
natural_cor_entropy_RSA_values = []
natural_intercept_values = []
natural_slope_values = []

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
all_temp_slope_values = []

temps = [0.0, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2]
modified_temps = [-0.2, 0.03, 0.1, 0.3, 0.6, 0.9, 1.2, 1.4]
temp_array = array(temps)
modified_temp_array = array(modified_temps)
protein_file = open("graph_mean_data_ordered_natural_noah.csv", "r")
natural_protein_data = protein_file.readlines()
protein_file.close()

header = natural_protein_data.pop(0)
for line in natural_protein_data:
    data = re.split("\t", line)
    natural_data.append(data)

for data in natural_data:
    #print data
    pdb_names.append(data[0])
    chain_names.append(data[1])
    natural_mean_RSA_values.append(data[2])
    natural_mean_entropy_values.append(data[3])
    natural_cor_entropy_RSA_values.append(data[4])
    natural_intercept_values.append(data[6])
    natural_slope_values.append(data[7])

natural_mean_RSA_values_array = analysis_functions.make_array(natural_mean_RSA_values)
natural_mean_entropy_values_array = analysis_functions.make_array(natural_mean_entropy_values)
natural_cor_entropy_RSA_values_array = analysis_functions.make_array(natural_cor_entropy_RSA_values)
natural_intercept_values_array = analysis_functions.make_array(natural_intercept_values)
natural_slope_values_array = analysis_functions.make_array(natural_slope_values)

protein_file_name = "graph_mean_data_ordered_0.0_noah.csv"
[designed_mean_RSA_values_00, designed_mean_entropy_values_00, designed_cor_entropy_RSA_values_00, designed_mean_KL_values_00, designed_intercept_values_00, designed_slope_values_00] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_00 = array(designed_mean_RSA_values_00)
designed_mean_entropy_values_array_00 = array(designed_mean_entropy_values_00)
designed_cor_entropy_RSA_values_array_00 = array(designed_cor_entropy_RSA_values_00)  
designed_mean_KL_values_array_00 = array(designed_mean_KL_values_00)
designed_intercept_values_array_00 = array(designed_intercept_values_00)
designed_slope_values_array_00 = array(designed_slope_values_00)

protein_file_name = "graph_mean_data_ordered_0.3_noah.csv"
[designed_mean_RSA_values_03, designed_mean_entropy_values_03, designed_cor_entropy_RSA_values_03, designed_mean_KL_values_03, designed_intercept_values_03, designed_slope_values_03] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_03 = array(designed_mean_RSA_values_03)
designed_mean_entropy_values_array_03 = array(designed_mean_entropy_values_03)
designed_cor_entropy_RSA_values_array_03 = array(designed_cor_entropy_RSA_values_03) 
designed_mean_KL_values_array_03 = array(designed_mean_KL_values_03)
designed_intercept_values_array_03 = array(designed_intercept_values_03)
designed_slope_values_array_03 = array(designed_slope_values_03)


protein_file_name = "graph_mean_data_ordered_0.6_noah.csv"
[designed_mean_RSA_values_06, designed_mean_entropy_values_06, designed_cor_entropy_RSA_values_06, designed_mean_KL_values_06, designed_intercept_values_06, designed_slope_values_06] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_06 = array(designed_mean_RSA_values_06)
designed_mean_entropy_values_array_06 = array(designed_mean_entropy_values_06)
designed_cor_entropy_RSA_values_array_06 = array(designed_cor_entropy_RSA_values_06)
designed_mean_KL_values_array_06 = array(designed_mean_KL_values_06) 
designed_intercept_values_array_06 = array(designed_intercept_values_06)
designed_slope_values_array_06 = array(designed_slope_values_06)

protein_file_name = "graph_mean_data_ordered_0.9_noah.csv"
[designed_mean_RSA_values_09, designed_mean_entropy_values_09, designed_cor_entropy_RSA_values_09, designed_mean_KL_values_09, designed_intercept_values_09, designed_slope_values_09] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_09 = array(designed_mean_RSA_values_09)
designed_mean_entropy_values_array_09 = array(designed_mean_entropy_values_09)
designed_cor_entropy_RSA_values_array_09 = array(designed_cor_entropy_RSA_values_09) 
designed_mean_KL_values_array_09 = array(designed_mean_KL_values_09)
designed_intercept_values_array_09 = array(designed_intercept_values_09)
designed_slope_values_array_09 = array(designed_slope_values_09)

protein_file_name = "graph_mean_data_ordered_1.2_noah.csv"
[designed_mean_RSA_values_12, designed_mean_entropy_values_12, designed_cor_entropy_RSA_values_12, designed_mean_KL_values_12, designed_intercept_values_12, designed_slope_values_12] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_12 = array(designed_mean_RSA_values_12)
designed_mean_entropy_values_array_12 = array(designed_mean_entropy_values_12)
designed_cor_entropy_RSA_values_array_12 = array(designed_cor_entropy_RSA_values_12) 
designed_mean_KL_values_array_12 = array(designed_mean_KL_values_12)
designed_intercept_values_array_12 = array(designed_intercept_values_12)
designed_slope_values_array_12 = array(designed_slope_values_12)

protein_file_name = "graph_mean_data_ordered_1.8_noah.csv"
[designed_mean_RSA_values_003, designed_mean_entropy_values_18, designed_cor_entropy_RSA_values_18, designed_mean_KL_values_18, designed_intercept_values_18, designed_slope_values_18] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_18 = array(designed_mean_RSA_values_18)
designed_mean_entropy_values_array_18 = array(designed_mean_entropy_values_18)
designed_cor_entropy_RSA_values_array_18 = array(designed_cor_entropy_RSA_values_18) 
designed_mean_KL_values_array_18 = array(designed_mean_KL_values_18)
designed_intercept_values_array_18 = array(designed_intercept_values_18)
designed_slope_values_array_18 = array(designed_slope_values_18)

protein_file_name = "graph_mean_data_ordered_2.4_noah.csv"
[designed_mean_RSA_values_24, designed_mean_entropy_values_24, designed_cor_entropy_RSA_values_24, designed_mean_KL_values_24, designed_intercept_values_24, designed_slope_values_24] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_24 = array(designed_mean_RSA_values_24)
designed_mean_entropy_values_array_24 = array(designed_mean_entropy_values_24)
designed_cor_entropy_RSA_values_array_24 = array(designed_cor_entropy_RSA_values_24)  
designed_mean_KL_values_array_24 = array(designed_mean_KL_values_24)
designed_intercept_values_array_24 = array(designed_intercept_values_24)
designed_slope_values_array_24 = array(designed_slope_values_24)

protein_file_name = "graph_mean_data_ordered_soft_noah.csv"
[designed_mean_RSA_values_soft, designed_mean_entropy_values_soft, designed_cor_entropy_RSA_values_soft, designed_mean_KL_values_soft, designed_intercept_values_soft, designed_slope_values_soft] = analysis_functions.get_mean_ordered_designed_data(protein_file_name)

designed_mean_RSA_values_array_soft = array(designed_mean_RSA_values_soft)
designed_mean_entropy_values_array_soft = array(designed_mean_entropy_values_soft)
designed_cor_entropy_RSA_values_array_soft = array(designed_cor_entropy_RSA_values_soft)  
designed_mean_KL_values_array_soft = array(designed_mean_KL_values_soft)
designed_intercept_values_array_soft = array(designed_intercept_values_soft)
designed_slope_values_array_soft = array(designed_slope_values_soft)

all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_00)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_03)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_06)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_09)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_12)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_18)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_24)
all_temp_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_soft)
all_temp_cor_entropy_RSA_values.append(natural_cor_entropy_RSA_values)


all_temp_slope_values.append(designed_slope_values_00)
all_temp_slope_values.append(designed_slope_values_03)
all_temp_slope_values.append(designed_slope_values_06)
all_temp_slope_values.append(designed_slope_values_09)
all_temp_slope_values.append(designed_slope_values_12)
all_temp_slope_values.append(designed_slope_values_18)
all_temp_slope_values.append(designed_slope_values_24)
all_temp_slope_values.append(designed_slope_values_soft)
all_temp_slope_values.append(natural_slope_values)

all_temp_cor_entropy_RSA_values_length = len(all_temp_cor_entropy_RSA_values[0])
all_temp_cor_entropy_RSA_values_array = []
all_temp_slope_values_array = []
index = 0
rcParams['font.size'] = 20
rcParams['figure.figsize'] = [14,6]
rcParams['lines.markersize'] = 8
fig = plt.figure(count, dpi = 500)
ax1 = axes([0.075,0.11,0.42, 0.85])
ax1.text(0.15, 1.0, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
#grid()
p1 = plot(designed_mean_entropy_values_array_00, designed_slope_values_array_00, color = "#EEEC8B", linestyle = "None", marker = "o") #Medium Yellow
p2 = plot(designed_mean_entropy_values_array_03, designed_slope_values_array_03, color = "#F8FC05", linestyle = "None", marker = "o") #Dark Yellow
p3 = plot(designed_mean_entropy_values_array_06, designed_slope_values_array_06, color = "#F5A112", linestyle = "None", marker = "o") #Dark Yellow/Dark Orange
p4 = plot(designed_mean_entropy_values_array_09, designed_slope_values_array_09, color = "#DF7401", linestyle = "None", marker = "o") #Orange
p5 = plot(designed_mean_entropy_values_array_12, designed_slope_values_array_12, color = "#D04900", linestyle = "None", marker = "o") #Dark Orange/ Dark Red
p6 = plot(designed_mean_entropy_values_array_18, designed_slope_values_array_18, color = "#FE2E2E", linestyle = "None", marker = "o") #Light Red
p7 = plot(designed_mean_entropy_values_array_24, designed_slope_values_array_24, color = "#B40404", linestyle = "None", marker = "o") #Dark Red
p8 = plot(designed_mean_entropy_values_array_soft, designed_slope_values_array_soft, color = "#800080", linestyle = "None", marker = "o") #Dark Red
p9 = plot(natural_mean_entropy_values_array, natural_slope_values_array, color = "#128EF5", linestyle = "None", marker = "o")
get_plot_format_entropy("Mean Entropy", "Slope", ax1)

ax2 = axes([0.57,0.11,0.42, 0.85])
text(-0.06, 1.0, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
p1 = plot(designed_intercept_values_array_00, designed_slope_values_array_00, color = "#EEEC8B", linestyle = "None", marker = "o") #Medium Yellow
p2 = plot(designed_intercept_values_array_03, designed_slope_values_array_03, color = "#F8FC05", linestyle = "None", marker = "o") #Dark Yellow
p3 = plot(designed_intercept_values_array_06, designed_slope_values_array_06, color = "#F5A112", linestyle = "None", marker = "o") #Dark Yellow/Dark Orange
p4 = plot(designed_intercept_values_array_09, designed_slope_values_array_09, color = "#DF7401", linestyle = "None", marker = "o") #Orange
p5 = plot(designed_intercept_values_array_12, designed_slope_values_array_12, color = "#D04900", linestyle = "None", marker = "o") #Dark Orange/ Dark Red
p6 = plot(designed_intercept_values_array_18, designed_slope_values_array_18, color = "#FE2E2E", linestyle = "None", marker = "o") #Light Red
p7 = plot(designed_intercept_values_array_24, designed_slope_values_array_24, color = "#B40404", linestyle = "None", marker = "o") #Dark Red
p8 = plot(designed_intercept_values_array_soft, designed_slope_values_array_soft, color = "#800080", linestyle = "None", marker = "o") #Dark Red
p9 = plot(natural_intercept_values_array, natural_slope_values_array, color = "#128EF5", linestyle = "None", marker = "o")
#This creates the figure legend
l1 = ax2.legend([p8, p1, p2, p3, p4, p5, p6, p7, p9],["Soft", "Fixed", "T = 0.3", "T = 0.6", "T = 0.9", "T = 1.2", "T = 1.8", "T = 2.4", "Natural"], "lower left", frameon = False, numpoints = 1, ncol = 2, prop = {'size': 15}) #Only one dot in the legend, no frame so "False", two columns so ncol = 2, and the font size is 15.

ax2.spines['right'].set_color('none')
#setp(ax2.get_yticklabels(), visible = False)
#ax2.yaxis.set_tick_params(size = 0)
get_plot_format_intercepts("Intercept", "Slope", ax2)

save_figure("Slope_Combination_Plot_Noah")
count = count + 1

#show()
