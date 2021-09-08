# -*- coding: utf-8 -*-
# MIT License

# Copyright (c) 2021 Dr. William A. Pisani

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Created on Sun Dec  1 16:10:15 2019

@author: wapisani

This code will plot pressure as a function of volume.
"""

# Import necessary libraries
import log
import os
import numpy as np
import matplotlib.pyplot as plt

# Define functions here
def moving_average(a, n):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


# Change to directory
directory = r"C:\Users\RDEL1WAP\Documents\Research\ERDC-EL-Git\Projects\Nylon6_CNT\Molecular_Dynamics\Composite"
os.chdir(directory)
filename = r"N6CNT1021IFFANPlyANEq08S3IE.log.lammps"
print(filename)
# Open log file
thermo = log.log(filename)

# Get pressure and volume data from log data
pressure, volume = thermo.get("Press","Volume")
# print(len(pressure))
pressure = pressure[1:]
volume = volume[1:]
# Compute average pressure and average volume during the two stages
total_time = 500000 # time of each stage, not total
thermo_freq = 10000 - 1
num_thermo_lines = int(total_time/thermo_freq)
press_avg_1 = np.average(pressure[:num_thermo_lines])
press_avg_2 = np.average(pressure[num_thermo_lines:])
vol_avg_1 = np.average(volume[:num_thermo_lines]) 
vol_avg_2 = np.average(volume[num_thermo_lines:])
#press_avg_1 = np.average(pressure[:num_thermo_lines])
#press_avg_2 = np.average(pressure[num_thermo_lines:])
#vol_avg_1 = np.average(volume[:num_thermo_lines])
#vol_avg_2 = np.average(volume[num_thermo_lines:])

# Compute bulk modulus
bulk_mod = -1*((press_avg_2-press_avg_1)/((vol_avg_2-vol_avg_1)/vol_avg_1)) *101325*10**-9 # in GPa

# Plot pressure-volume curve
fig, ax = plt.subplots()
ax.scatter(volume[1:],pressure[1:],marker='o',facecolor='none',edgecolors='#244395')
ax.scatter([vol_avg_1,vol_avg_2],[press_avg_1,press_avg_2],marker='^',facecolor='red',label="Avg press/vol")
ax.set_xlabel('Volume, Angstrom$^{3}$')
ax.set_ylabel('Pressure, atm')
ax.annotate(str(np.round(press_avg_1,1))+" atm",xy=(vol_avg_1-0.04*vol_avg_1,press_avg_1+1000))
ax.annotate(str(np.round(press_avg_2,1))+" atm",xy=(vol_avg_2+0.005*vol_avg_2,press_avg_2+500))
ax.legend()
ax.set_title("Bulk modulus = " + str(np.round(bulk_mod,4)) + " GPa")
fig.savefig(filename.split('.log')[0]+"_P_V_Plot.pdf")

print(f"Average low pressure: {press_avg_1}")
print(f"Average high pressure: {press_avg_2}")

# Try to access interaction energy
if filename.find("GNP") > -1:
    
    try:
        poly_gnp, gnp1_gnp2 = thermo.get("f_poly_gnp","f_gnp1_gnp2_avg")
        print(f"IE btwn polymer and matrix: {poly_gnp[-1]:3.3f}")
        print(f"IE btwn GNP1 and GNP2: {gnp1_gnp2[-1]:3.3f}")
    except:
        pass    