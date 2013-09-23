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
    frequency_data = []
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
    if count == 0: #Need to find a way to exclude the point
        frequency_data = [-1] # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #return frequency_data
    else:
        frequencies = [np.mean(bin_1)/np.mean(bin_1), np.mean(bin_2)/np.mean(bin_1), np.mean(bin_3)/np.mean(bin_1), np.mean(bin_4)/np.mean(bin_1), np.mean(bin_5)/np.mean(bin_1), np.mean(bin_6)/np.mean(bin_1), np.mean(bin_7)/np.mean(bin_1), np.mean(bin_8)/np.mean(bin_1), np.mean(bin_9)/np.mean(bin_1), np.mean(bin_10)/np.mean(bin_1), np.mean(bin_11)/np.mean(bin_1), np.mean(bin_12)/np.mean(bin_1), np.mean(bin_13)/np.mean(bin_1), np.mean(bin_14)/np.mean(bin_1), np.mean(bin_15)/np.mean(bin_1), np.mean(bin_16)/np.mean(bin_1), np.mean(bin_17)/np.mean(bin_1), np.mean(bin_18)/np.mean(bin_1), np.mean(bin_19)/np.mean(bin_1), np.mean(bin_20)/np.mean(bin_1)]
        for element in frequencies:
            if element == 0.0:
                frequency_data.append(0.0)
            else:
                frequency_data.append(np.log(element))
                #print np.log(element)
    #if (mean(bin_1)) == 0.0:
    #    print "MEAN OF BIN 1 is ZERO!!!!"
    #print frequency_data
    #frequency_data = [mean(bin_1), mean(bin_2), mean(bin_3), mean(bin_4), mean(bin_5)]
    #print "Number of residues in bin: " + str(count)
    return frequency_data

k = []
r = []
def set_k(RSA_midpoints):
    global k
    RSA_length = len(RSA_midpoints)
    for i in xrange(0,RSA_length):
        for element in k_start:
            k.append(element)
    k = np.array(k)

def set_r(RSA_midpoints):
    global r
    RSA_length = len(RSA_midpoints)
    for value in RSA_midpoints:
        for i in xrange(0,20):
            #print i
            r.append(value)
    r = np.array(r)

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
    #print k
    #m,k = get_x_(k_start)
    #a, b = params
    a = vars[0]
    b = vars[1]
    return (- (a + b*r)*k) - y_n

def lambda_function(exp_lambda1,y1, k):
    lam = exp_lambda1
    return (-lam*k) - y1


file = "align_natural_data_array_ordered_1QMV_A.dat"
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
all_midpoints = [((0.0 + 0.1)/2.0), ((0.1 + 0.2)/2.0), ((0.2 + 0.3)/2.0), ((0.3 + 0.4)/2.0), ((0.4 + 0.5)/2.0), ((0.5 + 0.6)/2.0), ((0.6 + 0.7)/2.0), ((0.7 + 0.8)/2.0), ((0.8 + 0.9)/2.0), ((0.9 + 1.0)/2.0)]
reg_midpoints =  [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
all_frequencies_list = [freq_bin_1, freq_bin_2, freq_bin_3, freq_bin_4, freq_bin_5, freq_bin_6, freq_bin_7, freq_bin_8, freq_bin_9, freq_bin_10]
midpoint_count = 0

for bin in all_frequencies_list:
    if(len(bin) != 1):
        freqs.append(bin)
        RSA_midpoints.append(all_midpoints[midpoint_count])
        RSA_points.append(reg_midpoints[midpoint_count])
    midpoint_count = midpoint_count + 1

all_frequencies = np.array(freqs)
set_r(RSA_midpoints)
set_k(RSA_midpoints)
x = get_x(RSA_midpoints)
#y_n = np.array(get_yn(log_freqs))
y_n = np.array(get_yn(all_frequencies))
#print k
#print r
#print y_n

vars = [1.09, -1.02]
#popt, pcov = curve_fit(function, x, y_n)
#popt = leastsq(residual, vars, args = (y_n, k, r))
#print popt
exp_lambda1 = 1.15
y1 = y_n[0:20]
print freq_bin_1
print np.exp(freq_bin_1)
popt_1 = leastsq(lambda_function, exp_lambda1, args = (freq_bin_1, k_start))
lambda1 = popt_1[0]
print lambda1

exp_lambda2 = 1
y2 = y_n[20:40]
popt_2 = leastsq(lambda_function, exp_lambda2, args = (freq_bin_2, k_start))
lambda2 = popt_2[0]
print lambda2

exp_lambda3 = 0.6
y3 = y_n[40:60]
popt_3 = leastsq(lambda_function, exp_lambda3, args = (freq_bin_3, k_start))
lambda3 = popt_3[0]
print lambda3

exp_lambda4 = 0.6
#y3 = y_n[40:60]
popt_4 = leastsq(lambda_function, exp_lambda4, args = (freq_bin_4, k_start))
lambda4 = popt_4[0]
print lambda4

exp_lambda5 = 0.6
#y3 = y_n[40:60]
popt_5 = leastsq(lambda_function, exp_lambda5, args = (freq_bin_5, k_start))
lambda5 = popt_5[0]
print lambda5

exp_lambda6 = 0.6
#y3 = y_n[40:60]
popt_6 = leastsq(lambda_function, exp_lambda6, args = (freq_bin_6, k_start))
lambda6 = popt_6[0]
print lambda6

exp_lambda7 = 0.5
#y3 = y_n[40:60]
popt_7 = leastsq(lambda_function, exp_lambda7, args = (freq_bin_7, k_start))
lambda7 = popt_7[0]
print lambda7

exp_lambda8 = 0.4
#y3 = y_n[40:60]
popt_8 = leastsq(lambda_function, exp_lambda8, args = (freq_bin_8, k_start))
lambda8 = popt_8[0]
print lambda8

exp_lambda9 = 0.3
#y3 = y_n[40:60]
popt_9 = leastsq(lambda_function, exp_lambda9, args = (freq_bin_9, k_start))
lambda9 = popt_9[0]
print lambda9

exp_lambda10 = 0.2
#y3 = y_n[40:60]
popt_10 = leastsq(lambda_function, exp_lambda10, args = (freq_bin_10, k_start))
lambda10 = popt_10[0]
print lambda10


#[slope, intercept] = popt
fig1 = figure(1)
x1 = linspace(0,18,num = 20)
y1 = - lambda1*x1*k_start
xlabel("k")
ylabel("log(frequency)")
#p1 = plt.plot(k_start,freq_bin_1, marker = 'o', color = 'blue', linestyle = 'None')
p2 = plt.plot(k_start,np.exp(freq_bin_1), marker = 'o', color = 'blue', linestyle = 'None')
#p3 = plt.plot(x1,y1, marker = 'o', color = 'blue', linestyle = '-')

fig2 = figure(2)
p1 = plt.plot(RSA_midpoints[0],lambda1, marker = 'o', color = 'black', linestyle = 'None')
p2 = plt.plot(RSA_midpoints[1],lambda2, marker = 'o', color = 'black', linestyle = 'None')
p3 = plt.plot(RSA_midpoints[2],lambda3, marker = 'o', color = 'black', linestyle = 'None')
p4 = plt.plot(RSA_midpoints[3],lambda4, marker = 'o', color = 'black', linestyle = 'None')
p5 = plt.plot(RSA_midpoints[4],lambda5, marker = 'o', color = 'black', linestyle = 'None')
p6 = plt.plot(RSA_midpoints[5],lambda6, marker = 'o', color = 'black', linestyle = 'None')
p7 = plt.plot(RSA_midpoints[6],lambda7, marker = 'o', color = 'black', linestyle = 'None')
p8 = plt.plot(RSA_midpoints[7],lambda8, marker = 'o', color = 'black', linestyle = 'None')
p9 = plt.plot(RSA_midpoints[8],lambda9, marker = 'o', color = 'black', linestyle = 'None')
p10 = plt.plot(RSA_midpoints[9],lambda10, marker = 'o', color = 'black', linestyle = 'None')
ylabel("Lambda")
xlabel("RSA")
plt.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]) 
plt.show()
