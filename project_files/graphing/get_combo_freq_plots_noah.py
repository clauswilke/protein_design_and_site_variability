import sys, os, math, string 
import analysis_functions as af
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

#Open all of the data files
#Get the frequencies
#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

#hydro_list = ['ILE', 'VAL', 'LEU', 'PHE', 'CYS', 'MET', 'ALA', 'GLY', 'THR', 'SER', 'TRP', 'TYR', 'PRO', 'HIS', 'GLU', 'GLN', 'ASP', 'ASN', 'LYS', 'ARG']
hydro_list = ['I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'S', 'W', 'Y', 'P', 'H', 'E', 'Q', 'D', 'N', 'K', 'R']
temps = ["0.0", "0.3", "0.6", "0.9", "1.2", "1.8", "2.4", "Soft"]
AA = resdict.values()

#This is dictionary that maps the natural sequence identifier with its PDB for the Noah Dataset
identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076':'2X1B', \
	'PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37', \
	'PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY', 'PF00355' :'1FQT', \
	'PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS', 'PF00542' :'1CTF', \
	'PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8',  \
	'PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI', 'PF02823' :'1AQT', \
	'PF04002':'2QLC','PF07686' :'2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

PDBS = identity_dict.values()
def get_mean_AA_freqs(file):
    freqs = []
    aa_counts = []
    input = open(file, "r")
    protein_data = input.readlines()
    input.close()
    
    site_data = af.get_transformed_data(protein_data)
    for site in site_data:
        new_site = []
        for element in site:
          new_site.append(float(element))
        freqs.append(new_site)
    freq_array = np.array(freqs)
 
    m,n = freq_array.shape
    total_aa = sum(freq_array)
    m,n = freq_array.shape

    j = 0
    while (j < n):
        aa_freq_sum = sum(freq_array[:, j])
        aa_counts.append(float(aa_freq_sum)/float(total_aa))
        j = j + 1
                       
    return aa_counts

def get_position_count_data(file):
    buried_counts = []
    surface_counts = []
    buried_sites = []
    surface_sites= []
    input = open(file, "r")
    protein_data = input.readlines()
    input.close()

    RSA = af.get_RSA_Values(file)
    alignment_length = len(RSA)
    site_data = af.get_transformed_data(protein_data)
    #print RSA
    #print site_data
    i = 0
    for site in site_data:
        if(float(RSA[i]) <= 0.05):
            buried_sites.append(site)
        else:
            surface_sites.append(site)
        i = i + 1

    buried_array = np.array(buried_sites)
    surface_array = np.array(surface_sites)

    buried_m, buried_n = buried_array.shape
    buried_total_sum = sum(sum(buried_array))
    surface_m, surface_n = surface_array.shape
    surface_total_sum = sum(sum(surface_array))
    #print surface_total_sum

    j = 0
    while (j < buried_n):
        buried_site_sum = sum(buried_array[:, j])
        buried_counts.append(float(buried_site_sum)/float(buried_total_sum))
        j = j + 1
    
    j = 0
    while (j < surface_n):
        surface_site_sum = sum(surface_array[:, j])
        surface_counts.append(float(surface_site_sum)/float(surface_total_sum))
        j = j + 1
    
    return buried_counts, surface_counts

rcParams['font.size'] = 14
rcParams['figure.figsize'] = [14,10]
#rcParams['lines.markersize'] = 8

N = 20 #The number of bars in each group - should be 20 for the 20 amino acids
index = np.arange(N) # The x locations for the groups 

count = 0

for temp in temps:
	designed_total_aa = []
	natural_total_aa = []
	designed_total_buried_aa = []
	designed_total_surface_aa = []
	natural_total_surface_aa = []
	natural_total_buried_aa = []
	for pdb in PDBS:
		pdb_id = pdb + "_A"
		designed_file = "align_data_array_" + pdb_id + "_" + str(temp) + ".dat" 
		natural_file = "align_natural_data_array_" + pdb_id  + ".dat" 
		#print designed_file
		designed_aa_freqs = get_mean_AA_freqs(designed_file)
		natural_aa_freqs = get_mean_AA_freqs(natural_file)
		designed_total_aa.append(designed_aa_freqs)
		natural_total_aa.append(natural_aa_freqs)
		
		[designed_buried_aa, designed_surface_aa] = get_position_count_data(designed_file)
		designed_total_buried_aa.append(designed_buried_aa)
		designed_total_surface_aa.append(designed_surface_aa)
		[natural_buried_aa, natural_surface_aa] = get_position_count_data(natural_file)
		natural_total_buried_aa.append(natural_buried_aa)
		natural_total_surface_aa.append(natural_surface_aa)
		
	designed_buried_array = np.array(designed_total_buried_aa)
	designed_surface_array = np.array(designed_total_surface_aa)
	designed_b_rows, designed_b_cols = designed_buried_array.shape
	designed_s_rows, designed_s_cols = designed_surface_array.shape
	designed_buried_mean_freqs = []
	designed_surface_mean_freqs = []
	
	natural_buried_array = np.array(natural_total_buried_aa)
	natural_surface_array = np.array(natural_total_surface_aa)
	natural_b_rows, natural_b_cols = natural_buried_array.shape
	natural_s_rows, natural_s_cols = natural_surface_array.shape
	natural_buried_mean_freqs = []
	natural_surface_mean_freqs = []
	
	designed_aa_freq_array = np.array(designed_total_aa)
	designed_aa_rows, designed_aa_cols = designed_aa_freq_array.shape
	designed_mean_aa_freqs = []
	natural_aa_freq_array = np.array(natural_total_aa)
	natural_aa_rows, natural_aa_cols = natural_aa_freq_array.shape
	natural_mean_aa_freqs = []

	for i in xrange(0, designed_b_cols):
		designed_buried_mean_freqs.append(np.mean(designed_buried_array[:, i]))
	for i in xrange(0, designed_s_cols):
		designed_surface_mean_freqs.append(np.mean(designed_surface_array[:,i]))
		
	for i in xrange(0, natural_b_cols):
		natural_buried_mean_freqs.append(np.mean(natural_buried_array[:, i]))
	for i in xrange(0, natural_s_cols):
		natural_surface_mean_freqs.append(np.mean(natural_surface_array[:,i]))			
	for i in xrange(0, designed_aa_cols):
		designed_mean_aa_freqs.append(np.mean(designed_aa_freq_array[:, i]))
	for i in xrange(0, natural_aa_cols):
		natural_mean_aa_freqs.append(np.mean(natural_aa_freq_array[:, i])) 
    
	if temp == "0.0":
		temp_string = "Designed, FB"
	else:
		temp_string = "Designed (Temp = " + str(temp) + ")"
		
	
	designed_buried_mean_freqs_dict = dict(zip(AA, designed_buried_mean_freqs))
	designed_surface_mean_freqs_dict = dict(zip(AA, designed_surface_mean_freqs))
	
	natural_buried_mean_freqs_dict = dict(zip(AA, natural_buried_mean_freqs))
	natural_surface_mean_freqs_dict = dict(zip(AA, natural_surface_mean_freqs))
		 
	designed_mean_aa_freqs_dict = dict(zip(AA, designed_mean_aa_freqs))	
	natural_mean_aa_freqs_dict = dict(zip(AA, natural_mean_aa_freqs))
	
	new_natural_buried_mean_freqs = []
	new_natural_surface_mean_freqs = []
	new_designed_buried_mean_freqs = []
	new_designed_surface_mean_freqs = []
	new_designed_mean_aa_freqs = []
	new_natural_mean_aa_freqs = []
	
	amino_acid_list = []
	
	for type in hydro_list:
		#amino = resdict[val
		#amino_acid_list.append(resdict[type])
		new_natural_buried_mean_freqs.append(natural_buried_mean_freqs_dict[type])
		new_natural_surface_mean_freqs.append(natural_surface_mean_freqs_dict[type])
		
		new_designed_buried_mean_freqs.append(designed_buried_mean_freqs_dict[type])
		new_designed_surface_mean_freqs.append(designed_surface_mean_freqs_dict[type])		
		
		new_designed_mean_aa_freqs.append(designed_mean_aa_freqs_dict[type])
		new_natural_mean_aa_freqs.append(natural_mean_aa_freqs_dict[type])

	#print temp_string
	#print "Designed Surface", designed_surface_mean_freqs
	#print "Natural Surface", natural_surface_mean_freqs
	fig = plt.figure(count)#, dpi = 500) 
	
	ax = plt.axes([0.09, 0.71, 0.87, 0.28])
	ax.text(0.32, 0.23, "All", ha = 'center', fontsize = 15)#, va = 'center', fontsize = 16)
	ax.text(-1.0, 0.251, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 15)
	width = 0.40 #The width of the bars
	b1 = plt.bar(index, new_designed_mean_aa_freqs, width, color = "blue")
	b2 = plt.bar(index+width, new_natural_mean_aa_freqs, width, color = "red")
	
	ax.set_xticklabels(hydro_list)
	ax.set_xticks(index+width)
	#ax.set_xlabel("Amino Acid")
	ax.set_ylabel("Frequency")
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25])
	ax.set_yticklabels(["0.00", "0.05", "0.10", "0.15", "0.20", "0.25"])
	ax.legend([b1, b2], [temp_string, "Natural"], numpoints = 1, frameon = False, loc = 1, prop = {'size': 15})
	
	ax2 = plt.axes([0.09, 0.38, 0.87, 0.28])
	ax2.text(0.8, 0.23, "Exposed", ha = 'center', fontsize = 15)
	ax2.text(-1.0, 0.251, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 15)
	b3 = plt.bar(index, new_designed_surface_mean_freqs, width, color = "blue")
	b4 = plt.bar(index + width, new_natural_surface_mean_freqs, width, color = "red") 
	ax2.set_xticklabels(hydro_list)
	ax2.set_xticks(index+width)
	#ax2.set_xlabel("Amino Acid")
	ax2.set_ylabel("Frequency")
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax2.get_xaxis().tick_bottom()
	ax2.get_yaxis().tick_left()
	ax2.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25])
	ax2.set_yticklabels(["0.00", "0.05", "0.10", "0.15", "0.20", "0.25"])
	
	ax3 = plt.axes([0.09, 0.06, 0.87, 0.28])
	ax3.text(0.64, 0.23, "Buried", ha = 'center', fontsize = 15)
	ax3.text(-1.0, 0.251, "C", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 15)
	b5 = plt.bar(index, new_designed_buried_mean_freqs, width, color = "blue")
	b6 = plt.bar(index + width, new_natural_buried_mean_freqs, width, color = "red") 
	
	ax3.set_xticklabels(hydro_list)
	ax3.set_xticks(index+width)
	ax3.set_xlabel("Amino Acid")
	ax3.set_ylabel("Frequency")
	ax3.spines['top'].set_visible(False)
	ax3.spines['right'].set_visible(False)
	ax3.get_xaxis().tick_bottom()
	ax3.get_yaxis().tick_left()
	ax3.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25])
	ax3.set_yticklabels(["0.00", "0.05", "0.10", "0.15", "0.20", "0.25"])
	#plt.legend([b5, b6], [temp_string, "Natural"], numpoints = 1, frameon = False, loc = 1 , prop = {'size': 15})
	
	fig_title = "Noah_Freq_Combo_Plots_" + str(temp) + ".pdf"
	plt.savefig(fig_title, format = None)
	#plt.show()
	count = count + 1
