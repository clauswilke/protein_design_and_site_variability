import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
import numpy as np
import analysis_functions as af
from scipy.optimize import curve_fit
#import matplotlib as plt
from pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import leastsq 
import matplotlib.cm as cm
import random 
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
    #natural_distribution = af.get_AA_distribution(natural_proteins)
    natural_distribution = get_AA_distribution_mod(natural_proteins)
    #natural_dis_array = array(natural_distribution)
    #m,n = natural_dis_array.shape
    #print "num_residues, length of alignment: " + str(n),m
    #print natural_distribution
    natural_RSA = af.get_RSA_Values(natural_proteins)
    natural_RSA_array = af.make_array(natural_RSA)
    seq_length = len(natural_RSA)
    frequency_data = []
    all_k_values = []
    k_values = []
    k_lists = []
    bin_1 = []
    bin_2 = []
    bin_3 = []
    bin_4 = []
    bin_5 = []
    bin_6 = []
    bin_7 = []
    bin_8 = []
    bin_9 = []
    bin_10 = []
    bin_11 = []
    bin_12 = []
    bin_13 = []
    bin_14 = []
    bin_15 = []
    bin_16 = []
    bin_17 = []
    bin_18 = []
    bin_19 = []
    bin_20 = []
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
            bin_6.append(site[5])
            bin_7.append(site[6])
            bin_8.append(site[7])
            bin_9.append(site[8])
            bin_10.append(site[9])
            bin_11.append(site[10])
            bin_12.append(site[11])
            bin_13.append(site[12])
            bin_14.append(site[13])
            bin_15.append(site[14])
            bin_16.append(site[15])
            bin_17.append(site[16])
            bin_18.append(site[17])
            bin_19.append(site[18])
            bin_20.append(site[19])
            i = i + 1
            count = count + 1
        else:
            i = i + 1
    if count != 0: #Need to find a way to exclude the point
        frequencies = [np.mean(bin_1)/np.mean(bin_1), np.mean(bin_2)/np.mean(bin_1), np.mean(bin_3)/np.mean(bin_1), np.mean(bin_4)/np.mean(bin_1), np.mean(bin_5)/np.mean(bin_1), np.mean(bin_6)/np.mean(bin_1), np.mean(bin_7)/np.mean(bin_1), np.mean(bin_8)/np.mean(bin_1), np.mean(bin_9)/np.mean(bin_1), np.mean(bin_10)/np.mean(bin_1), np.mean(bin_11)/np.mean(bin_1), np.mean(bin_12)/np.mean(bin_1), np.mean(bin_13)/np.mean(bin_1), np.mean(bin_14)/np.mean(bin_1), np.mean(bin_15)/np.mean(bin_1), np.mean(bin_16)/np.mean(bin_1), np.mean(bin_17)/np.mean(bin_1), np.mean(bin_18)/np.mean(bin_1), np.mean(bin_19)/np.mean(bin_1), np.mean(bin_20)/np.mean(bin_1)]
        all_k_values = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
        k_count = 0
        for element in frequencies:
            if element == 0.0:
                continue
            else:
                frequency_data.append(np.log(element))
                k_values.append(all_k_values[k_count])
                k_count = k_count + 1
    set_k_r(k_values, lower_RSA_boundary, upper_RSA_boundary)
    return frequency_data


def set_k(k_values):
    global k
    for value in k_values:
        k.append(k_values)
    
def set_k_r(k_values, lower_RSA_boundary, upper_RSA_boundary):
    global k
    global r
    global k_lists
    for value in k_values:
        k.append(value)
        r.append((lower_RSA_boundary + upper_RSA_boundary)/2.0)
    k_lists.append(k_values)

def get_x(RSA_midpoints):
    #global k
    r = []
    RSA_length = len(RSA_midpoints)
    for value in RSA_midpoints:
        for i in xrange(0,20):
            #print i
            r.append(value)
    x = np.array(r)*np.array(k)
    return x

def get_yn(log_freqs):
    y_n = []
    for value in log_freqs:
        for freq in value:
            y_n.append(freq)
    return y_n

def residual(vars, y_n, k, r):
    a = vars[0]
    b = vars[1]
    d = (-(a + b*r)*k) - y_n
#    print d
#    print "RMSD: ", np.sqrt(sum(d*d))
    #print k
    #print r

    return d

def myleastsq(fun, vars, args):
    res = fun( vars, args[0], args[1], args[2])
    d = sum(res*res)
    boolean = True
    while boolean == True:
        new_vars = [vars[0] + random.gauss(0, 0.1), vars[1] +  random.gauss(0, 0.1)]
        new_d = sum(fun( new_vars, args[0], args[1], args[2] )**2)
        #print "New Vars"
        #print new_vars
        #print np.sqrt(d), np.sqrt(new_d)
        if new_d < d:
            vars=new_vars
            if abs( d-new_d )/d < .0001:
                boolean = False
            d = new_d
            #break
    #print vars
    return vars


#This is dictionary that maps the natural sequence identifier with its PDB for the Noah Dataset
identity_dict = {'PF00013' :'1WVN','PF00018' : '2O9S','PF00041':'1TEN','PF00072' :'1MVO','PF00076':'2X1B', \
	'PF00085' :'1FB0','PF00111' :'1CZP','PF00168' :'3F04','PF00169' :'1UNQ','PF00179' :'1Z2U','PF00226' :'2O37', \
	'PF00240' :'2BWF','PF00249' :'1GUU','PF00254' :'2PPN','PF00313' :'3I2Z','PF00327' :'1BXY', 'PF00355' :'1FQT', \
	'PF00364' :'2EVB','PF00381' :'1PTF','PF00439' :'3JVL','PF00486' :'2ZXJ','PF00498' :'3GQS', 'PF00542' :'1CTF', \
	'PF00550' :'1T8K','PF00581' :'1GN0','PF00582' :'2Z3V','PF00595' :'2H3L','PF00691' :'1OAP','PF00708' :'3BR8',  \
	'PF01029' :'1TZV','PF01035' :'3GVA','PF01451' :'1JL3','PF01627' :'2A0B','PF01833' :'3MQI', 'PF02823' :'1AQT', \
	'PF04002':'2QLC','PF07686' :'2PND','PF08666' :'1UCS','PF12844' :'3FYM','PF07679' :'1U2H'}

PDBS = identity_dict.values()
pp = PdfPages("Noah_PDB_Linear_Model_Plots.pdf")
count = 0
chain = "A"
for pdb_id in PDBS:
	file = "align_natural_data_array_ordered_" + pdb_id + "_" + chain_id + ".dat"
	k = []
	r = []
	k_lists = []
	print file
	#file = "align_natural_data_array_ordered_1B4T_A.dat"
	freq_bin_1 = get_RSA_frequencies(file, 0.0, 0.1)
	freq_bin_2 = get_RSA_frequencies(file, 0.1, 0.2) 
	freq_bin_3 = get_RSA_frequencies(file, 0.2, 0.3)
	freq_bin_4 = get_RSA_frequencies(file, 0.3, 0.4) 
	freq_bin_5 = get_RSA_frequencies(file, 0.4, 0.5)
	freq_bin_6 = get_RSA_frequencies(file, 0.5, 0.6) 
	freq_bin_7 = get_RSA_frequencies(file, 0.6, 0.7) 
	freq_bin_8 = get_RSA_frequencies(file, 0.7, 0.8) 
	freq_bin_9 = get_RSA_frequencies(file, 0.8, 0.9) 
	freq_bin_10 = get_RSA_frequencies(file, 0.9, 1.0)

	k_start = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
	freqs = [] 
	RSA_midpoints = []
	RSA_points = []
	RSA_midpoints = [((0.0 + 0.1)/2.0), ((0.1 + 0.2)/2.0), ((0.2 + 0.3)/2.0), ((0.3 + 0.4)/2.0), ((0.4 + 0.5)/2.0), ((0.5 + 0.6)/2.0), ((0.6 + 0.7)/2.0), ((0.7 + 0.8)/2.0), ((0.8 + 0.9)/2.0), ((0.9 + 1.0)/2.0)]
	reg_midpoints =  [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
	all_frequencies_list = [freq_bin_1, freq_bin_2, freq_bin_3, freq_bin_4, freq_bin_5, freq_bin_6, freq_bin_7, freq_bin_8, freq_bin_9, freq_bin_10]
	all_frequencies = np.array(all_frequencies_list)
	k = np.array(k)
	r = np.array(r)
	y_n = np.array(get_yn(all_frequencies))
	vars = [1.09, -1.02]
	#popt, pcov = curve_fit(function, x, y_n)
	#popt = leastsq(residual, vars, args = (y_n, k, r))
	solution = myleastsq(residual, vars, args = (y_n, k, r))
	rcParams['font.size'] = 16
	intercept = solution[0]
	slope  = solution[1]
	print intercept
	print slope
	fig = figure(count)
	ax = plt.axes([0.12, 0.10, 0.85, 0.85]) 
	p1 = plot(k_lists[0], freq_bin_1, marker= 'o', linestyle = 'None', color = cm.hot(.1))
	p2 = plot(k_lists[1], freq_bin_2, marker= 'o', linestyle = 'None', color = cm.hot(.2))
	p3 = plot(k_lists[2], freq_bin_3, marker= 'o', linestyle = 'None', color = cm.hot(.3))
	p4 = plot(k_lists[3], freq_bin_4, marker= 'o', linestyle = 'None', color = cm.hot(.4))
	p5 = plot(k_lists[4], freq_bin_5, marker= 'o', linestyle = 'None', color = cm.hot(.5))
	p6 = plot(k_lists[5], freq_bin_6, marker= 'o', linestyle = 'None', color = cm.hot(.6))
	p7 = plot(k_lists[6], freq_bin_7, marker= 'o', linestyle = 'None', color = cm.hot(.7))
	p8 = plot(k_lists[7], freq_bin_8, marker= 'o', linestyle = 'None', color = cm.hot(.8))
	p9 = plot(k_lists[8], freq_bin_9, marker= 'o', linestyle = 'None', color = cm.hot(.9))
	p10 = plot(k_lists[9], freq_bin_10, marker= 'o', linestyle = 'None', color = cm.hot(1))


	knew = arange(0,20,.1)
	r1 = RSA_midpoints[0]
	r2 = RSA_midpoints[1]
	r3 = RSA_midpoints[2]
	r4 = RSA_midpoints[3]
	r5 = RSA_midpoints[4]
	r6 = RSA_midpoints[5]
	r7 = RSA_midpoints[6]
	r8 = RSA_midpoints[7]
	r9 = RSA_midpoints[8]
	r10 = RSA_midpoints[9]
	a = intercept
	b = slope

	ynew1 = -(a + b*r1)*knew
	ynew2 = -(a + b*r2)*knew
	ynew3 = -(a + b*r3)*knew
	ynew4 = -(a + b*r4)*knew
	ynew5 = -(a + b*r5)*knew
	ynew6 = -(a + b*r6)*knew
	ynew7 = -(a + b*r7)*knew
	ynew8 = -(a + b*r8)*knew
	ynew9 = -(a + b*r9)*knew
	ynew10 = -(a + b*r10)*knew

	p1 = plot(knew, ynew1, linestyle = '-', color = cm.hot(.1))
	p2 = plot(knew, ynew2, linestyle = '-', color = cm.hot(.2))
	p3 = plot(knew, ynew3, linestyle = '-', color = cm.hot(.3))
	p4 = plot(knew, ynew4, linestyle = '-', color = cm.hot(.4))
	p5 = plot(knew, ynew5, linestyle = '-', color = cm.hot(.5))
	p6 = plot(knew, ynew6, linestyle = '-', color = cm.hot(.6))
	p7 = plot(knew, ynew7, linestyle = '-', color = cm.hot(.7))
	p8 = plot(knew, ynew8, linestyle = '-', color = cm.hot(.8))
	p9 = plot(knew, ynew9, linestyle = '-', color = cm.hot(.9))
	p10 = plot(knew, ynew10, linestyle = '-', color = cm.hot(1))

	ylabel ("log(Frequency)")
	xlabel("Amino Acid Rank k")
	xlim(-0.25, 20)
	ylim(-12.25, 0)
	plt.xticks([0, 5, 10, 15, 20] )
	plt.yticks([-12, -10, -8, -6, -4, -2, 0])
	name = pdb_id + "_" + chain_id 
	slope_string = "Slope: %.3f" % (slope)
	intercept_string = "Intercept: %.3f" % (intercept)
	ax.text(0.25, -11, name, fontsize = 16)
	ax.text(0.25, -11.5, slope_string, fontsize = 16)
	ax.text(0.25, -12, intercept_string, fontsize = 16)
	pp.savefig(fig)
	count = count + 1
	#plt.show()

pp.close()
