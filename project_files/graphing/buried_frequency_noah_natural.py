import math, os, string, re
import analysis_functions as af
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

#This is dictionary that maps the natural sequence identifier with its PDB for the Noah Dataset
identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076':'2X1B', \
	'PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37', \
	'PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY', 'PF00355' :'1FQT', \
	'PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS', 'PF00542' :'1CTF', \
	'PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8',  \
	'PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI', 'PF02823' :'1AQT', \
	'PF04002':'2QLC','PF07686' :'2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

temps = ["0.0", "0.3","0.6", "0.9", "1.2", "1.8", "2.4"]
chain_id = 'A'
PDBS = identity_dict.values()

mean_buried_percentage = []
#The inputs to this functions are three lists with the RSA, entropy, and KL site data for a protein. This sorts the entropy, and 
#KL data by RSA into three categories: buried, partially buried, and surface. It then returns the data. 
def get_percent_buried(RSA_data):
    num_sites = len(RSA_data)
    i = 0
    num_buried = 0
    while i < num_sites:
        if (float(RSA_data[i])<0.05):
            num_buried = num_buried + 1
            i = i + 1
        else:
            i = i + 1
    percent_buried = (float(num_buried)/float (num_sites))
    return percent_buried

buried_percentage_list = []
#Open a pdf
#pp = PdfPages("Duncan_PDB_Linear_Frequency_Plots.pdf")
rcParams['font.size'] = 12
rcParams['figure.figsize'] = [14,6]
#rcParams['lines.markersize'] = 8

chain = 'A'
#pp = PdfPages("Duncan_Natural_Buried_Bar_Graphs.pdf")
count = 0

buried_percentage_list = []
for pdb_id in PDBS:
    file = "align_natural_data_array_" + pdb_id + "_" +  chain_id + ".dat"
    print "Processing file: " + file	
    natural_proteins = file
    natural_RSA = af.get_RSA_Values(natural_proteins)
    buried_percentage = get_percent_buried(natural_RSA)
    buried_percentage_list.append(buried_percentage)


N = 40 #The number of bars in each group - should be 20 for the 20 amino acids
index = np.arange(N) # The x locations for the groups 
fig = plt.figure(count, dpi = 500) 
ax = plt.axes([0.07, 0.14, 0.92, 0.83])
width = 0.45 #The width of the bars
b1 = plt.bar(index, buried_percentage_list, width, color = "blue")
ax.set_xticklabels(PDBS, rotation = 'vertical')
width_tick = (float (width)/2.0)
ax.set_xticks(index + width_tick)
ax.set_xlabel("PDB")
ax.set_ylabel("Frequency of Buried Residues")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.set_yticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45])
ax.set_yticklabels(["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45"])
fig_title = "Noah_Percent_Buried_Natural.pdf"
plt.savefig(fig_title, format = None)
#plt.show()


