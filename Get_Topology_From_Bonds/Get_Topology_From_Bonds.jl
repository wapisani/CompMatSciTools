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

# This program will read in the bonds from the designated file and then write out the 
# topology (angles, dihedrals, impropers) to a file.
# Future improvements
# Enable reading from a standard LAMMPS data file
# Enable writing bonds, angles, dihedrals, and impropers to a LAMMPS data file
# Checking if correct number of arguments are supplied

# Include necessary libraries and functions
using StatsBase
include("angles.jl")
include("dihedrals.jl")
include("impropers.jl")
include("write_topology.jl")

bond_dir = pwd()
example_file = ARGS[1]

function parsebonds(directory::String,file::String)
    open(joinpath(directory,file),"r") do f
        lines = split(read(f,String),"\n")
        bonds = Array{Int64}[]
        for line in lines
            line = split(line," ")
            if length(line) == 4 # If bond type is present
                atom1 = parse(Int64,line[3])
                atom2 = parse(Int64,line[4])
            elseif length(line) == 3 # If bond type is not present
                atom1 = parse(Int64,line[2])
                atom2 = parse(Int64,line[3])
            end
            push!(bonds,[atom1,atom2])
        end
        return bonds
    end
end

bonds = parsebonds(bond_dir,example_file)

# Get angles from bonds
angle_array = angles(bonds)

# Get dihedrals
dihedral_array = dihedrals(bonds,angle_array)

# Get impropers
improper_array = impropers(bonds)

outfile = join([splitext(example_file)[1],"_topology_test.dat"])
write_topology(bond_dir,outfile,bonds,angle_array,dihedral_array,improper_array)
