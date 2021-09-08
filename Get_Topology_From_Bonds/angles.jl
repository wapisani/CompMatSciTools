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
    angles(bonds::Array{Array{Int64,N} where N,1}) -> Array{Tuple{Int64,Int64,Int64},1}
    
Determines the angles of a molecule or molecular system from the array of bonds given as input. 
"""
function angles(bonds::Array{Array{Int64,N} where N,1})  # returns a set of tuples
    println("Determining angles... ")
    angle_set = Set(Tuple{Int64,Int64,Int64}[])
    count = 1
    for i in 1:length(bonds)
        atom1 = bonds[i][1]
        atom2 = bonds[i][2]
        atom12set = Set([atom1,atom2])
        for j in 1:length(bonds)
            atom3 = bonds[j][1]
            atom4 = bonds[j][2]
            atom34set = Set([atom3,atom4])
            atom1234set = Set([atom1,atom2,atom3,atom4])
            if length(atom1234set) == 3
                centeratom = collect(setdiff(atom12set,symdiff(atom12set,atom34set)))
                outsideatoms = collect(symdiff(atom12set,atom34set))
                angle = (outsideatoms[1],centeratom[1],outsideatoms[2])
                if ∉(angle,angle_set)
                    angle = reverse(angle)
                    push!(angle_set,angle)
                    count += 1
                end
            end
        end
    end
    println("Angles determined!")
    return collect(angle_set) # Return Array{Tuple{Int64,Int64,Int64},1}
end