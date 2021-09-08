#!/usr/bin/env python3
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
Created on Tue Sep  7 23:36:31 2021

@author: Dr. William A. Pisani

ElasticMD is a Python package I have been working on that will
help streamline the process of analyzing stress-strain data
generated from uniaxial tensile and shearing simulations in 
LAMMPS. This script (when finished) will show how ElasticMD
can be used to extract data from log files, export it to CSV
for analysis in third-party programs, compute elastic moduli 
from stress-strain data, and make publication-worthy figures.
"""

import os
import ElasticMD

os.chdir(r'./Examples/PEEK')
example_log_file = r'PEEK0517RxYM3CrS1_3.log.lammps'

# Extract thermodynamic data from example log file
thermo = ElasticMD.extract.log(example_log_file)

# Log file was from uniaxial tensile strain simulation
# in the z-direction so get the strain and stress in z-direction
eengz,szz_ave = thermo.get(('eengz','szz_ave'),0)





