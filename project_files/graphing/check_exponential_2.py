import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
from scipy.stats import *
import matplotlib.pyplot as plt
import analysis_functions as af
from matplotlib.backends.backend_pdf import PdfPages

#This function takes a list of sequences and returns a list of lists that contain all of the frequency data for every site
def get_transformed_data(list_of_sequences):
    transformed_distribution = []
    RSA = []
    for seq in list_of_sequences:
        sequences = seq.strip()
    list_of_sequences.pop(0)

    new_data = []
    for line in list_of_sequences:
        element = line.split()
        element.pop(0)
        element.pop(0)
        element.pop(0)
        new_data.append(element)
    for data in new_data:
        new_elements = []
        for count in data:       
            new_count = int(count)
            new_elements.append(new_count)
        transformed_distribution.append(new_elements)
    return transformed_distribution

def get_AA_distribution_mod(proteins):
    #print proteins
    protein_distribution = []
    transformed_protein_distribution = []
    num_AA = 0
    input = open(proteins, "r")
    protein_data = input.readlines()
    input.close()
    protein_distribution = get_transformed_data(protein_data)
    for site in protein_distribution:
        new_site = site
        num_AA = sum(new_site)
        aa_probs = []
        for count in new_site:   
            if count == 0:
                prob = 0.0
                #prob = float(1)/float(num_AA + 20)
            else:
                prob = float(count)/float(num_AA)
                #prob = float(count + 1)/float(num_AA + 20)
            aa_probs.append(prob)
        transformed_protein_distribution.append(aa_probs)
    return transformed_protein_distribution


def get_RSA_frequencies(natural_proteins, lower_RSA_boundary, upper_RSA_boundary):
    natural_distribution = af.get_AA_distribution(natural_proteins)
    #natural_distribution = get_AA_distribution_mod(natural_proteins)
    #natural_dis_array = array(natural_distribution)
    #m,n = natural_dis_array.shape
    #print "num_residues, length of alignment: " + str(n),m
    #print natural_distribution
    natural_RSA = af.get_RSA_Values(natural_proteins)
    natural_RSA_array = af.make_array(natural_RSA)
    seq_length = len(natural_RSA)
    bin_1 = []
    bin_2 = []
    bin_3 = []
    bin_4 = []
    bin_5 = []
    i = 0
    count = 0
    for site in natural_distribution:
        if (lower_RSA_boundary<=natural_RSA_array[i] and natural_RSA_array[i]<= upper_RSA_boundary):
            #print natural_RSA_array[i]
            #print site[0:4]
            bin_1.append(site[0])
            bin_2.append(site[1])
            bin_3.append(site[2])
            bin_4.append(site[3])
            bin_5.append(site[4])
            i = i + 1
            count = count + 1
        else:
            i = i + 1
    if count == 0:
        frequency_data = [0.0, 0.0, 0.0, 0.0, 0.0]
    else:
        frequency_data = [mean(bin_1)/mean(bin_1), mean(bin_2)/mean(bin_1), mean(bin_3)/mean(bin_1), mean(bin_4)/mean(bin_1), mean(bin_5)/mean(bin_1)]
    if (mean(bin_1)) == 0.0:
        print "MEAN OF BIN 1 is ZERO!!!!"
    print frequency_data
    #frequency_data = [mean(bin_1), mean(bin_2), mean(bin_3), mean(bin_4), mean(bin_5)]
    print "Number of residues in bin: " + str(count)
    return frequency_data

file = "align_natural_data_array_ordered_2BR9_A.dat"
#natural_freq_array_1 = array(get_RSA_frequencies(file, 0.0, 0.1))
#natural_freq_array_2 = array(get_RSA_frequencies(file, 0.1, 0.2))
#natural_freq_array_3 = array(get_RSA_frequencies(file, 0.2, 0.3))
#natural_freq_array_4 = array(get_RSA_frequencies(file, 0.3, 0.4))
#natural_freq_array_5 = array(get_RSA_frequencies(file, 0.4, 0.5))
#natural_freq_array_6 = array(get_RSA_frequencies(file, 0.5, 0.6))
#natural_freq_array_7 = array(get_RSA_frequencies(file, 0.6, 0.7))
#natural_freq_array_8 = array(get_RSA_frequencies(file, 0.7, 0.8))
#natural_freq_array_9 = array(get_RSA_frequencies(file, 0.8, 0.9))
#natural_freq_array_10 = array(get_RSA_frequencies(file, 0.9, 1.0))

#print natural_freq_array_1
#print natural_freq_array_2
#print natural_freq_array_3
#print natural_freq_array_4
#print natural_freq_array_5
#print natural_freq_array_6
#print natural_freq_array_7
#print natural_freq_array_8
#print natural_freq_array_9
#print natural_freq_array_10


#This searches for all of the unaligned sequences 
searchStr = "align_natural_data_array_ordered_" + "[a-zA-Z0-9_\.\-]*" +  ".dat"
pp = PdfPages("Duncan_PDB_Linear_Frequency_Plots.pdf")
count = 0
for path, names, filename in os.walk('.',False): #Searchs for all the *.dat files for the #natural protein 
    for file in filename:
        if(re.search(searchStr, file)!=None):
            fileparts = re.split("_",file)
            pdb_id = fileparts[5].upper() #Gets the PDB Name and the chain_id
            chain_id = fileparts[6]
            chain_id = chain_id[0]
            
            natural_proteins = file
            print natural_proteins

            #Calculates all of the data for comparison (ex. entropy)
            natural_freq_array_1 = array(get_RSA_frequencies(file, 0.0, 0.1))
            natural_freq_array_2 = array(get_RSA_frequencies(file, 0.2, 0.3))
            natural_freq_array_3 = array(get_RSA_frequencies(file, 0.5, 0.6))
            natural_freq_array_4 = array(get_RSA_frequencies(file, 0.8, 0.9))
            #print natural_freq_array_1
            rcParams['font.size'] = 16
            rcParams['figure.figsize'] = [10, 6]
            rcParams['lines.linewidth'] = 2.0
            rcParams['lines.markersize'] = 5
            rcParams['lines.marker'] = 'o'            
            fig = plt.figure(count, dpi = 400)
            ax = axes([0.12, 0.10, 0.80, 0.85])
            p1 = ax.plot([1,2,3,4,5],log(natural_freq_array_1), color = 'black', linestyle = '-', marker = 'o')
            p2 = ax.plot([1,2,3,4,5],log(natural_freq_array_2), color = 'blue', linestyle = '-', marker = 'o')
            p3 = ax.plot([1,2,3,4,5],log(natural_freq_array_3), color = 'red', linestyle = '-', marker = 'o')
            p4 = ax.plot([1,2,3,4,5], log(natural_freq_array_4), color = 'green', linestyle = '-', marker = 'o')
            #This creates the figure legend
            text(5.25, 0, pdb_id, ha = 'center', va = 'center', fontsize = 16)
            l1 = ax.legend([p1, p2, p3, p4],['RSA = [0.0, 0.1]','RSA = [0.2, 0.3]' ,'RSA = [0.5, 0.6]','RSA = [0.8, 0.9]'], "lower left", frameon = False, numpoints = 1, ncol = 1, prop = {'size': 15}) #Only one dot in the legend, no frame so "False", two columns so ncol = 2, and the font size is 15.
            xlabel("Amino Acid")
            ylabel("Relative Frequencies")
            plt.yticks()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            plt.yticks([log(0.003), log(0.01), log(0.02), log(0.05), log(0.1), log(0.2), log(0.5), log(1.0)], ["0.003", "0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1.0"])
            plt.xticks([1, 2, 3, 4, 5], ["1", "2", "3", "4", "5"])
            plt.xlim(0.9,5.1)
            #plt.ylim(-4.61,0.0)
            plt.ylim(log(0.002),0.0)
            pp.savefig(fig)
            count = count + 1
           
pp.close()

