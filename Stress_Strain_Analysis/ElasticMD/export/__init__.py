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
Created on Wed Nov 25 16:30:25 2020

@author: RDEL1WAP
"""

def export2Csv(lists2Export,csvname,csvdir,header=''):
    """
    Takes a list of lists and exports them to a CSV file called "csvname"
    located in csvdir.

    Returns
    -------
    None.

    """
    import os
    import numpy as np
    os.chdir(csvdir)
    ncolumns = len(lists2Export) 
    nrows = len(lists2Export[0])
    
    # Add lists to NumPy 2D array
    data = np.zeros((nrows,ncolumns))
    for col in range(ncolumns):
        data[:,col] = lists2Export[col]
    
    # Save NumPy data
    np.savetxt(csvname,data,delimiter=",",header=header)
    print(f"{csvname} saved to {csvdir}")