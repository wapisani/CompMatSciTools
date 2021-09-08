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
Created on Tue Sep  1 14:39:10 2020

@author: William A. Pisani, Ph.D. 

This module will read in a LAMMPS log file and extract all thermodynamic data.
"""

import numpy as np

class log:
    
    def __init__(self,logname,suppressOutput=False):
        self.keywords = [] # List of lists of thermo keywords (Step, Temp, Press, etc.), if different runs have different numbers of keywords
        self.nkeywords = [] # List of numbers of keywords
        self.keywordIndices = [] # List of dictionaries of thermo keyword-index pairs corresponding to each run
        self.headerIndices = [] # List of indices where thermo keywords are
        self.endRunIndices = [] # List of indices where "Loop time string is, marking the end of each run
        self.data = [] # List of 2d numpy arrays
        self.logname = logname # Name of log file
        self.extractKeywordsAndData()
        
        if suppressOutput == False:
            self.printHeaders()
    
    def __repr__(self):
        return f"{self.__class__.__name__} object containing thermodynamic data from {self.logname}"
    
    def printHeaders(self):
        """
        Prints thermodynamic headers in order.
        
        Returns
        -------
        None.

        """
        print(f"Data extraction from {self.logname} was successful!")
        print(f"{len(self.keywords)} thermodynamic section(s) found")
        print("\nThermodynamic keyword headers are as follows (index, header):")
        print("-------------------------------------------------------------")
        for index,header in enumerate(self.keywords):
            print(f"{index}\t{header}\n")
        if len(self.keywords) > 1:
            print("Please note that since more than one (1) section of non-identical thermodynamic data was found, you will need to specify which section of data you wish to extract.")
            print("For example, thermo.get(('Step','Temp','Press'),0) to get the step, temperature, and pressure from the first section of data")
    
    
    def extractKeywordsAndData(self):
        """
        Get all keyword/data sections from all LAMMPS runs in log file 

        Returns
        -------
        None.

        """
        with open(self.logname,'r') as logfile:
            logContents = logfile.read()
            
        splitLogContents = logContents.split('\n')
        
        for index,line in enumerate(splitLogContents):
            if line.find("Per MPI ") > -1 or line.find("Memory ") > -1: # Per MPI is for modern versions of LAMMPS, Memory is for older versions
                self.headerIndices.append(index+1) # Per MPI rank memory always occurs one line before the thermo keywords line
            elif line.find("Loop time ") > -1:
                self.endRunIndices.append(index)
        
        for index,headerIndex in enumerate(self.headerIndices):
            # Thermo keywords
            line = splitLogContents[headerIndex]
            headerLine = " ".join(line.split()).split(' ')
            
            # Raw Data
            start, stop = headerIndex+1, self.endRunIndices[index]
            rawData = splitLogContents[start:stop]
            
            if headerLine not in self.keywords:
                self.keywords.append(headerLine)
                self.nkeywords.append(len(headerLine))
                
                keywordPairs = {}
                for index,keyword in enumerate(headerLine):
                    keywordPairs.update({keyword:index})
                self.keywordIndices.append(keywordPairs)
                
                # Convert raw data to numpy array
                npData = np.zeros((len(rawData),len(headerLine)))
                for i,dataLine in enumerate(rawData):
                    dataLine = " ".join(dataLine.split())
                    for j,value in enumerate(dataLine.split(' ')):
                        npData[i,j] = value
                self.data.append(npData)
                
            else: # If thermo header is identical to one already stored
                # Get index of first occurence of thermo header that is identical to the next header/data set to be stored
                
                firstOccurenceIndex = self.keywords.index(headerLine)
                # Convert raw data to numpy array
                npData = np.zeros((len(rawData),len(headerLine)))
                for i,dataLine in enumerate(rawData):
                    dataLine = " ".join(dataLine.split())
                    for j,value in enumerate(dataLine.split(' ')):
                        npData[i,j] = value
                # Add numpy array to numpy array of first occurence
                self.data[firstOccurenceIndex] = np.concatenate((self.data[firstOccurenceIndex],npData))
            
 
    def get(self,keys,index=0):
        """
        Parameters
        ----------
        keys : tuple
            Tuple of thermodynamic keywords. Common examples are "Step", "Temp", and "Press".
            Fixes, computes, and variables can also be extracted if given the appropriate term (e.g. "f_sxx_ave").
        index : int, optional
            Index for which section of data you wish to pull from. The default is 0.

        Raises
        ------
        Exception
            If no keywords are specified, an exception will be raised. At least one keyword must be specified.

        Returns
        -------
        properties : list
            List of 1d numpy arrays corresponding to the input tuple.

        """
        if len(keys) == 0:
            raise Exception("no keywords specified, you must specify at least one keyword (e.g. Step, Temp, etc)")
        
        if type(keys) == tuple:
            properties = []
            for key in keys:
                keyValue = self.keywordIndices[index][key]
                data = self.data[index][:,keyValue]
                properties.append(data)
        elif type(keys) == str:
            keyValue = self.keywordIndices[index][keys]
            data = self.data[index][:,keyValue]
            properties = data
        
        return properties
            


    