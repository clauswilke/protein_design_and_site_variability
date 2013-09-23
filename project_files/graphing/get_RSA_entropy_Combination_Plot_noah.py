import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
from scipy.stats import *
import matplotlib.pyplot as plt
import analysis_functions
from matplotlib.backends.backend_pdf import PdfPages

#This function takes a list of string values and then returns a numpy array of float values
def make_array(list_of_values):
    new_value_list = []
    for value in list_of_values:
        new_value = value.rstrip()
        #print new_value
        new_value = float(new_value)
        new_value_list.append(new_value)
    value_array = array(new_value_list)
    return value_array

def format_RSA_entropy_combination_plot():
    plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0], ["0.0", "1.0", "2.0", "3.0", "4.0"])
    plt.xticks([0, 0.25, 0.5, 0.75, 1], ["0", "", "0.5", "", "1"])
    plt.xlim(-0.075,1.0)
    plt.ylim(-0.25,4.0)


#This is dictionary that maps the natural sequence identifier with its PDB for the Noah Dataset
identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076':'2X1B', \
	'PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37', \
	'PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY', 'PF00355' :'1FQT', \
	'PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS', 'PF00542' :'1CTF', \
	'PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8',  \
	'PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI', 'PF02823' :'1AQT', \
	'PF04002':'2QLC','PF07686' :'2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

pp = PdfPages("Noah_PDB_Combination_Plots.pdf")
PDBS = identity_dict.values()
count = 1 #Count used to plot figures
for protein in PDBS:
    chain_id = 'A'
    pdb_id = protein
    #pdb_part = protein[0:4].upper()
    #chain_part = protein[4]
    file = "graph_data_" + pdb_id + "_" + chain_id + "_natural_noah.csv"
    natural_proteins = file #These lines get the chain, pdb_id filename, 
    #fileparts = re.split("_", file)
    #pdb_id = fileparts[2].upper()
    #chain_id = fileparts[3]

    natural_file = open(file)
    natural_file_data = natural_file.readlines()
    natural_data = []
    for line in natural_file_data:
        data = re.split("\t", line)
        natural_data.append(data)

    natural_entropy = natural_data[0]
    natural_entropy.pop(0)
    natural_entropy_array = analysis_functions.make_array(natural_entropy)  

    natural_RSA = natural_data[1]
    natural_RSA.pop(0)
    natural_RSA_array = analysis_functions.make_array(natural_RSA)

    designed_filename_00 = "graph_data_" + pdb_id + "_" + chain_id + "_0.0_noah.csv"
    [designed_entropy_array_00, designed_RSA_array_00, designed_KL_array_00] = analysis_functions.get_designed_graph_data(designed_filename_00)
    designed_filename_03 = "graph_data_" + pdb_id + "_" + chain_id + "_0.3_noah.csv"
    [designed_entropy_array_03, designed_RSA_array_03, designed_KL_array_03] = analysis_functions.get_designed_graph_data(designed_filename_03)
    designed_filename_06 = "graph_data_" + pdb_id + "_" + chain_id + "_0.6_noah.csv"
    [designed_entropy_array_06, designed_RSA_array_06, designed_KL_array_06] = analysis_functions.get_designed_graph_data(designed_filename_06)
    designed_filename_09 = "graph_data_" + pdb_id + "_" + chain_id + "_0.9_noah.csv"
    [designed_entropy_array_09, designed_RSA_array_09, designed_KL_array_09] = analysis_functions.get_designed_graph_data(designed_filename_09)
    designed_filename_12 = "graph_data_" + pdb_id + "_" + chain_id + "_1.2_noah.csv"
    [designed_entropy_array_12, designed_RSA_array_12, designed_KL_array_12] = analysis_functions.get_designed_graph_data(designed_filename_12)
    designed_filename_18 = "graph_data_" + pdb_id + "_" + chain_id + "_1.8_noah.csv"
    [designed_entropy_array_18, designed_RSA_array_18, designed_KL_array_18] = analysis_functions.get_designed_graph_data(designed_filename_18)
    designed_filename_24 = "graph_data_" + pdb_id + "_" + chain_id + "_2.4_noah.csv"
    [designed_entropy_array_24, designed_RSA_array_24, designed_KL_array_24] = analysis_functions.get_designed_graph_data(designed_filename_24)


    rcParams['font.size'] = 16
    rcParams['figure.figsize'] = [10, 6]
    rcParams['lines.linewidth'] = 2.0
    rcParams['lines.markersize'] = 5
    rcParams['lines.marker'] = 'o'

    #This figure plots Correlation of RSA, Entropy at Sites vs Mean Entropy
    fig = plt.figure(count, dpi = 400)
    ax = axes([0.08, 0.15, 0.20, 0.35])
    p1 = ax.plot(designed_RSA_array_12, designed_entropy_array_12, color = "black", linestyle = "None", marker = "o")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ylabel("Entropy")
    xlabel("RSA")
    format_RSA_entropy_combination_plot()
    text(0.17, 3.6, "T = 1.2", ha = 'center', va = 'center', fontsize = 16)           

    ax2 = axes([0.31, 0.15, 0.20, 0.35])
    p2 = ax2.plot(designed_RSA_array_18, designed_entropy_array_18, color = "black", linestyle = "None", marker = "o")
    setp(ax2.get_yticklabels(), visible = False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    xlabel("RSA")
    format_RSA_entropy_combination_plot()
    text(0.17, 3.6, "T = 1.8", ha = 'center', va = 'center', fontsize = 16)

    ax3 = axes([0.54, 0.15, 0.20, 0.35])
    p3 = ax3.plot(designed_RSA_array_24, designed_entropy_array_24, color = "black", linestyle = "None", marker = "o")
    setp(ax3.get_yticklabels(), visible = False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.get_xaxis().tick_bottom()
    ax3.get_yaxis().tick_left()
    xlabel("RSA")
    format_RSA_entropy_combination_plot()
    text(0.17, 3.6, "T = 2.4", ha = 'center', va = 'center', fontsize = 16)

    ax4 = axes([0.77, 0.15, 0.20, 0.35])
    p4 = ax4.plot(natural_RSA_array, natural_entropy_array, color = "black", linestyle = "None", marker = "o")
    setp(ax4.get_yticklabels(), visible = False)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.get_xaxis().tick_bottom()
    ax4.get_yaxis().tick_left()
    xlabel("RSA")
    format_RSA_entropy_combination_plot()
    text(0.17, 3.6, "Natural", ha = 'center', va = 'center', fontsize = 16)

    ax5 = axes([0.08, 0.55, 0.20, 0.35])
    p5 = ax5.plot(designed_RSA_array_00, designed_entropy_array_00, color = "black", linestyle = "None", marker = "o")
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    setp(ax5.get_xticklabels(), visible = False)
    ax5.get_xaxis().tick_bottom()
    ax5.get_yaxis().tick_left()
    ylabel("Entropy")
    format_RSA_entropy_combination_plot()
    text(0.11, 3.6, "Fixed", ha = 'center', va = 'center', fontsize = 16)

    ax6 = axes([0.31, 0.55, 0.20, 0.35])
    p6 = ax6.plot(designed_RSA_array_03, designed_entropy_array_03, color = "black", linestyle = "None", marker = "o")
    setp(ax6.get_yticklabels(), visible = False)
    setp(ax6.get_xticklabels(), visible = False)
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    ax6.get_xaxis().tick_bottom()
    ax6.get_yaxis().tick_left()
    format_RSA_entropy_combination_plot()
    text(0.205, 3.6, "T = 0.3", ha = 'center', va = 'center', fontsize = 16)

    ax7 = axes([0.54, 0.55,0.20, 0.35])
    p7 = ax7.plot(designed_RSA_array_06, designed_entropy_array_06, color = "black", linestyle = "None", marker = "o")
    setp(ax7.get_yticklabels(), visible = False)
    setp(ax7.get_xticklabels(), visible = False)
    ax7.spines['top'].set_visible(False)
    ax7.spines['right'].set_visible(False)
    ax7.get_xaxis().tick_bottom()
    ax7.get_yaxis().tick_left()
    format_RSA_entropy_combination_plot()
    text(0.17, 3.6, "T = 0.6", ha = 'center', va = 'center', fontsize = 16)

    ax8 = axes([0.77, 0.55, 0.20, 0.35])
    p8 = ax8.plot(designed_RSA_array_09, designed_entropy_array_09, color = "black", linestyle = "None", marker = "o")
    setp(ax8.get_yticklabels(), visible = False)
    setp(ax8.get_xticklabels(), visible = False)
    ax8.spines['top'].set_visible(False)
    ax8.spines['right'].set_visible(False)

    ax8.get_xaxis().tick_bottom()
    ax8.get_yaxis().tick_left()
    format_RSA_entropy_combination_plot()
    text(0.17, 3.6, "T = 0.9", ha = 'center', va = 'center', fontsize = 16)

    fig_title = "KL-Divergence vs Entropy"
    save_fig_title = "RSA_vs_Entropy_" + pdb_id + "_Combination_Plot_Noah" +  ".pdf"
    savefig(save_fig_title, format = None)
    pp.savefig(fig)
    count = count + 1
    print pdb_id
    #if(pdb_id == '1TEN'):
    #    show()
    #    quit()

pp.close()
