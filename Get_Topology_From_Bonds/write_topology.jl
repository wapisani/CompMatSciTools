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
    write_topology(directory::String,file::String,bonds::Array{Array{Int64,N} where N,1},angles::Array{Tuple{Int64,Int64,Int64},1},dihedrals::Array{NTuple{4,Int64},1},impropers::Array{NTuple{4,Int64},1})

Writes out the topology to a user-specified text file
"""
function write_topology(directory::String,file::String,bonds::Array{Array{Int64,N} where N,1},angles::Array{Tuple{Int64,Int64,Int64},1},dihedrals::Array{NTuple{4,Int64},1},impropers::Array{NTuple{4,Int64},1})
    outfile = open(joinpath(directory,file),"w")
    write(outfile,"# Please note that the types of the topology have not been determined and have not been included here.\n\n")
    nbonds = length(bonds)
    nangles = length(angles)
    ndihedrals = length(dihedrals)
    nimpropers = length(impropers)
    write(outfile,"$nbonds bonds\n$nangles angles\n$ndihedrals dihedrals\n$nimpropers impropers\n\n")
    count = 1
    write(outfile,"Bonds \n\n")
    for bond in bonds
        atom1 = bond[1]
        atom2 = bond[2]
        write(outfile,"$count $atom1 $atom2\n")
        count += 1
    end
    count = 1
    write(outfile,"\nAngles\n\n")
    for angle in angles
        atom1 = angle[1]
        atom2 = angle[2]
        atom3 = angle[3]
        write(outfile,"$count $atom1 $atom2 $atom3\n")
        count += 1
    end
    count = 1
    write(outfile,"\nDihedrals\n\n")
    for dihedral in dihedrals
        atom1 = dihedral[1]
        atom2 = dihedral[2]
        atom3 = dihedral[3]
        atom4 = dihedral[4]
        write(outfile,"$count $atom1 $atom2 $atom3 $atom4\n")
        count += 1
    end
    count = 1
    write(outfile,"\nImpropers\n\n")
    for improper in impropers
        atom1 = improper[1]
        atom2 = improper[2]
        atom3 = improper[3]
        atom4 = improper[4]
        write(outfile,"$count $atom1 $atom2 $atom3 $atom4\n")
        count += 1
    end
    close(outfile)
    println("Finished writing topology to $file in $directory")
end