#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created on Mon Feb 19 13:30:52 2018
#
#@author: ClayElmore
#"""

# This program will curve fit a section of a spectrum with a Gaussian 

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# -------------------- This is your section to mess with -------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# This is the file name. IT HAS TO BE EXACT.
#file_name = "PEGDA 1000 LiSPA 20 03_21_18 1.csv"
#file_name = "PEGDA 1000 01_02_18.csv"
#file_name = "PEGDA 1000 LiSS 20 03_19_18 1.csv"
#file_name = "PEGDA 1000 LiVS 20 03_19_18 1.csv"
#file_name = "PEGDA 1000 NaSTFSI 20 03_28_18 1.csv"
#file_name = "PEGDA 1000 KSTFSI 20 03_28_18 1.csv"
#file_name = "PEGDA 1000 LiTFSI 20 02_15_18 1.csv"
#file_name = "PEGDA 1000 MgTFSI 20 02_21 1.csv"
#file_name = "PEGDA 1500 NaTFSI 30 03_05_18 1.csv"
#file_name = "PEGDA 1500 LiTFSI 30 03_02_18 1.csv"
#file_name = "PEGDA 1500 KTFSI 30 03_01_18 1-1.csv"
#file_name = "PEGDA 1500 KTFSI 30 04_17_18 1.csv"
#file_name = "PEGDA 1500 CaSTFSI 30 04_17_18 1.csv"
file_name = "PEGDA 1500 MgTFSI 30 03_05_18 1.csv"
#plt.savefig('test.png',dpi=1000)

# This lets you analyze only the section of a spectra that you want
section_start = 718 # starting wavenumber
points = 47 # number of points analyzed 

# how many peaks do you think there are 
# This has to be less than or equal to your guesses
peak_number = 3

# wavenumber to scale on
scale_peak = 747
scale_amount = 1

# guesses for where you think the peaks are 
x_0_guess_1 = [747,734,729]

# guess for half peakwidths
lam_guess_1 = [5.9,4,4]

# guess for the height correction
y_guess = 0.2

# This will give the guesses for the peak intensities. 
# note: I am currently estimating these right now in "Section 4" of the code. 
# if you find that that they are just not working, comment out all of section
# 4 and just place the guesses here. if I find a better way to do section 4 I 
# will let you know
A_guess_1 = lam_guess_1
A_guess_2 = A_guess_1


# bounds for the equtions, these are quite important. GET THEM CLOSE 
bnds=[]
for i in range(peak_number):
    bnds.append((x_0_guess_1[i]-2,x_0_guess_1[i]+2))
    bnds.append((0.1,lam_guess_1[i]+10))
    bnds.append((A_guess_1[i]/50,A_guess_1[i]*100))
bnds.append((0,1))


# name for the final plot
title_name = file_name

# plot axes and legend names
x_axis = 'Wavenumbers (cm^-1)'
y_axis = 'Activities (Normalized)'
plot_legend = ['Experimental','Combined Peaks','Individual Peaks']
scale_legend = ['Pure','SS','TFSI']

# this section controls what you want to plot
full_plot = False # plots the original spectrum 
fit_plot = True # plots the Gaussian fitted curve
close_all = True # clears all previous plots
hold_on = False # plots everything on one for the scaled plots no fitting

# if the algorightm is not converging, try a different one in the list here
# Note: Some of these do not handle the bounds. I recommend mainly sticking 
# to SLSQP and L-BFGS-B but any of them should get you to the same result
choices = ['SLSQP','TNC','Nelder-Mead','Powell','CG','BFGS','L-BFGS-B','COBYLA']
algorithm = choices[0]

# if the algorithm needs more iternations you can increase it here
iterations = 1000

# Wavelet noise reduction:
wavelet_correction = False
wavelet_plot = False
wavelet_compare = True
epsilon = 0.01
level_trans = 5

calc_li = False

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# -------------------- PLEASE STOP HERE. STOOOOPPPPPPPPP -------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #


# DO NOT TOUCH THE REST OF THIS. PLEASE. I DON'T WANT TO MESS WITH THIS AGAIN.
# --------------------------------------------------------------------------- #
import math
import csv
from matplotlib import pyplot as plt
pi = math.pi
e = math.e
import numpy as np
from scipy.optimize import minimize
from pywt import wavedec as dwt
from pywt import waverec as dwt_inv

# --------------------------------------------------------------------------- #
# Section 1: Define the functions that will govern the curves

# Gaussian distribution for just one curve
def gaus_dist(x,x_01,lam1,A1,y):
    return (y+A1*(1/lam1/np.sqrt(2*pi)*np.exp(-(x-x_01)**2/2/(lam1**2))))

def gaus_dist_mult(x,x_01,lam1,A1,y):
    ans = 0
    for i in range(len(lam1)):
        ans += y+(A1[i]*(1/lam1[i]/np.sqrt(2*pi)*np.exp(-(x-x_01[i])**2/2/(lam1[i]**2)))) 
    return (ans)    

# --------------------------------------------------------------------------- #
# Section 2: Import the data 

def get_data(file_name):
    # initialize the arrays for the data
    wavenumbers = []
    activities = []
    
    # open the csv file
    with open(file_name) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        
        # These variables will allow you to only include the data
        row_count = 0;
        end = 0
        for row in readCSV:
            if row == []:
                end=1
            elif row_count > 18 and not end:
                activities.append(float(row[1]))
                wavenumbers.append(float(row[0]))
            row_count+=1
        return (wavenumbers,activities)

# run the above function
wavenumbers,activities = get_data(file_name)
wavenumbers_old = wavenumbers
activities_old = activities

# scale the data
def scale(activities,wavenumbers):
    wave_mod = []
    for i in range(len(wavenumbers)):
        wave_mod.append(math.floor(wavenumbers[i]))
        
    scale_index = wave_mod.index(scale_peak)
    scale_act = max(activities[scale_index-5:scale_index+5])
    for i in range(len(activities)):
        activities[i] = activities[i]/scale_act
    
    return(activities)
        
activities = scale(activities,wavenumbers)


# user controlled plots
if close_all:
    plt.close('all')
    
if wavelet_correction:
    activities_old = activities

    activities = np.array(activities)
    coeff = dwt(activities,'db4',level = level_trans)
    # make the wavelet:
    x_wavelet = []
    for i in coeff:
        for j in i:
            x_wavelet.append(j)
    
    # --------------------------------------------------------------------------- #
    # --------------------------------------------------------------------------- #
    # Noise Reduction
            
    def reduce_noise (decomp,epsilon):
        for i in range(len(decomp)):
            if i>0:
                for j in range(len(decomp[i])):
                    if abs(decomp[i][j]) < epsilon:
                        decomp[i][j] = 0
        return(decomp)
    
    coeff = reduce_noise(coeff,epsilon)
    
    # make the new vector
    activities = dwt_inv(coeff,'db4')
    
    # show the new wavelet vector
    x_wavelet_reduc = []
    for i in coeff:
        for j in i:
            x_wavelet_reduc.append(j)
            
    activities=activities[:-1]
    
    # --------------------------------------------------------------------------- #
    # --------------------------------------------------------------------------- #
    # Plotting
    
    activities = scale(activities,wavenumbers)
    
    if wavelet_plot:
        plt.figure()
        plt.plot(x_wavelet)
        plt.plot(x_wavelet_reduc,'--')
        plt.title('Wavelets')
        plt.legend(['Level '+str(level_trans)+' Decomposition',\
                    'Level '+str(level_trans)+' Decomposition W/ Noise Reduction'])
    
    if wavelet_compare:
        plt.figure()
        plt.plot(wavenumbers,activities_old)
        plt.plot(wavenumbers,activities)
        plt.legend(['Original','Reduced Noise'])
        plt.title('Wavelet Noise Reduction')

if full_plot and not hold_on:
    plt.plot(wavenumbers,activities,'k')
    plt.title(title_name)
    plt.xlabel('Wavenumber')
    plt.ylabel('Activities (A.U.)')
    
if full_plot and hold_on:
    plt.plot(wavenumbers,activities)
    plt.title(title_name)
    plt.xlabel('Wavenumber')
    plt.ylabel('Activities (A.U.)')
    plt.legend(scale_legend)
    
# --------------------------------------------------------------------------- #
# Section 3: Slice the data to whatever section you want 

def section(wavenumbers,activities,section_start,points): 
    
    # find the index of the place that the section starts
    wave_mod = []
    for i in range(len(wavenumbers)):
        wave_mod.append(math.floor(wavenumbers[i]))
        
    low_index = wave_mod.index(section_start)
    
    wavenumbers_mod = wavenumbers[low_index : low_index+points]
    activities_mod = activities[low_index : low_index+points]
    
    return (wavenumbers_mod,activities_mod)

# run the function
wavenumbers,activities = section(wavenumbers,activities,section_start,points)

# --------------------------------------------------------------------------- #
# Section 5: This section will curve fit the given section to a given number
#            of curves for however many peaks there are suspected to be

# make the x and y arrays out of the spectrum's defined section
x = np.array(wavenumbers)
y = np.array(activities)

# create the guess array as a numpy array
guess = []
for i in range(peak_number):
    guess.append(x_0_guess_1[i])
    guess.append(lam_guess_1[i])
    guess.append(A_guess_1[i])
guess.append(y_guess)
guess = np.array(guess)

# define the two norm of the Gaussian function:
def gaus_cauchy_dist_mult_two_norm(var,y,x):
    y_new=[]
    for j in range(len(y)):
        ans = 0
        for i in range(peak_number):
            ans += var[-1]+(var[i*3+2]*(1/var[i*3+1]/np.sqrt(2*pi)\
                    *np.exp(-(x-var[i*3])**2/2/(var[i*3+1]**2)))) 
        y_new.append(ans)
    
    y_new = np.array(y_new)
    res = np.linalg.norm(y-y_new)
    return (res)

# use the sciPy optimize function to fit the curve
bnds = tuple(bnds)
params = minimize(gaus_cauchy_dist_mult_two_norm, guess, args=(y, x),\
                  method=algorithm,bounds=bnds,options={'maxiter':iterations})

params =  params.x

# store the parameters
x_01 = []
lam1 = []
A1 = []
x_02 = []
lam2 = []
A2 = []
for i in range(peak_number):
    x_01.append(params[3*i])
    lam1.append(params[3*i+1])
    A1.append(params[3*i+2])
y_0 = (params[-1])
# --------------------------------------------------------------------------- #
# Section 7: Plotting

# plot the theoretical curve
if fit_plot:
    # plot the sliced data
    plt.figure()
    plt.plot(wavenumbers,activities,'mo')
    wave_est_new = np.linspace(section_start-10,section_start+10+points,500)
    act_est = []
    for i in wave_est_new:
        act_est.append(gaus_dist_mult(i,x_01,lam1,A1,y_0))
    plt.plot(wave_est_new,act_est,'g')
    
    
    # plot the individual peaks
    for i in range(peak_number):
        act_est1 = []
        for j in wave_est_new:
            act_est1.append(gaus_dist(j,x_01[i],lam1[i],A1[i],y_0))
        plt.plot(wave_est_new,act_est1,'k')
    
    plt.title(title_name,**{'size':'18','weight':'bold','fontname':'Arial'})
    plt.xlabel(x_axis,**{'size':'14','fontname':'Arial'})
    plt.ylabel(y_axis,**{'size':'14','fontname':'Arial'})
    plt.legend(plot_legend)

#def est_li(A_free,A_tethered):

def error_analysis (activities_old,activities,wavenumbers,A1,x_01,lam1,y):
#    # first find residuals:
#    fit_act = []
#    for i in wavenumbers:
#        fit_act.append(gaus_dist_mult(i,x_01,lam1,A1,y_0))
#    fit_act = np.array(fit_act)
#    activities = np.array(activities)
#    residuals = abs(activities - fit_act)
#    
#    # find the std and mean of the residuals 
#    residual_mean = np.mean(residuals)

    # take a section of spectrum that has no peaks and find a residual using this 
    noise_points = activities_old[len(activities_old)-100:len(activities_old)]
    noise_mean = np.mean(noise_points)
    residual = abs(noise_points - noise_mean)
    residual_mean = np.average(residual)
    print('Residual mean:',residual_mean)
    
    # monte carlo simulation:    
    # these are the vectors that hold all of the results
    num_trial = 1000
    A_test_tethered = []
    x0_test_tethered =[]
    lam_test_tethered = []
    A_test_free = []
    x0_test_free =[]
    lam_test_free = []
    y_test= []
    
    # this simulates all the data
    for i in range(num_trial):
        
        # create the random data
        act_mont = gaus_dist_mult(wavenumbers,x_01,lam1,A1,y)
        for j in range(len(act_mont)):
            act_mont[j] = act_mont[j] + residual_mean*np.random.normal()    
        
        # fit the simulated data
        coeffs = minimize(gaus_cauchy_dist_mult_two_norm, guess, args=(act_mont, wavenumbers),\
                  method=algorithm,bounds=bnds,options={'maxiter':iterations})
        coeffs =  coeffs.x
        
        # this stores the coefficients
#        x0_test_tethered.append(coeffs)
#        lam_test_tethered.append(coeffs[1])
        A_test_tethered.append(coeffs[2])            
#        x0_test_free.append(coeffs[3])
#        lam_test_free.append(coeffs[4])
        A_test_free.append(coeffs[5])
#        y_test.append(coeffs[-1])
    
    
    # calcualte the theoretical ion percentages for each of the tests
    cation_sim = []
    for i in range(len( A_test_free)):
        cation_sim.append(A_test_free[i] / (A_test_free[i]+A_test_tethered[i]))
    cation_error= np.std(cation_sim)*2
    free_cation_sim = np.mean(cation_sim)
    
    # calculate free ion percentage
    free_cation = A1[1] / (A1[1]+A1[0])
    
    return free_cation, free_cation_sim, cation_error

if calc_li:   
    free_cation, free_cation_sim, cation_error = \
        error_analysis (activities_old, activities,wavenumbers,A1,x_01,lam1,y_0)
    
    print('Experimental percent of free Li: '+str(free_cation*100))
    print('Simulation percent of free Li: '+str(free_cation_sim*100)\
          +' Plus or Minus '+str(cation_error*100))
    print('Percent Error: '+str(cation_error/free_cation*100))

if not calc_li:
    print('Free cation: ',(A1[1] / (A1[1]+A1[0]))*100)

