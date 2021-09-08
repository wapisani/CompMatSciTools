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
    dihedrals(bonds::Array{Array{Int64,N} where N,1},angles::Array{Tuple{Int64,Int64,Int64},1}) -> Array{NTuple{4,Int64},1}

Determines the dihedrals of a molecule or molecular system given the bonds and angles.
"""
function dihedrals(bonds::Array{Array{Int64,N} where N,1},angles::Array{Tuple{Int64,Int64,Int64},1})
    println("Determining dihedrals...")
    dihedral_set = Set(Tuple{Int64,Int64,Int64,Int64}[])
    count = 1
    for i = 1:length(angles)
        atom1 = angles[i][1]
        atom2 = angles[i][2]
        atom3 = angles[i][3]
        atom123set = Set([atom1,atom2,atom3])
        for j = 1:length(bonds)
            atom4 = bonds[j][1]
            atom5 = bonds[j][2]
            atom45set = Set([atom4,atom5])
            atom1345set = Set([atom1,atom3,atom4,atom5])
            atom12345set = Set([atom1,atom2,atom3,atom4,atom5])
            if length(atom12345set) == 4 && length(atom1345set) == 3
                commonatom = setdiff(atom45set,symdiff(atom123set,atom45set))
                newatom = collect(setdiff(atom45set,commonatom))
                if collect(commonatom)[1] == atom1
                    dihedral = (newatom[1],atom1,atom2,atom3)
                else
                    dihedral = (atom1,atom2,atom3,newatom[1])
                end
                if ∉(dihedral,dihedral_set)
                    dihedral = reverse(dihedral)
                    push!(dihedral_set,dihedral)
                    count += 1
                end
            end
        end
    end
    println("Dihedrals determined!")
    return collect(dihedral_set)
end