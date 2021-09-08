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
    impropers(bonds::Array{Array{Int64,N} where N,1}) -> Array{NTuple{4,Int64},1}

Determines the impropers of a molecule or molecular system given the bonds.
"""
function impropers(bonds::Array{Array{Int64,N} where N,1})
    println("Determining impropers...")
    count_atomid1s = countmap(atomid[1] for atomid in bonds)
    count_atomid2s = countmap(atomid[2] for atomid in bonds)
    count_atomids = mergewith(+,count_atomid1s,count_atomid2s)

    # Get atomids with 3 or more bonds, those atoms have impropers
    atomids = [key for (key,value) in count_atomids if value >= 3]

    improper_array = Tuple{Int64,Int64,Int64,Int64}[]

    # Find atom ids connected to improper atoms
    for atomid in atomids
        atomid_connections = []
        for bond in bonds
            if atomid == bond[1]
                push!(atomid_connections,bond[2])
            elseif atomid == bond[2]
                push!(atomid_connections,bond[1])
            end
        end

        # Now that connections are found, get combinations of impropers
        if length(atomid_connections) == 3 # If only three bonds
            imp = (atomid_connections[1],atomid,atomid_connections[2],atomid_connections[3])
            push!(improper_array,imp)

        elseif length(atomid_connections) == 4 # If four bonds
            imp = (atomid_connections[1],atomid,atomid_connections[2],atomid_connections[3]) # 1-2-3
            push!(improper_array,imp)

            imp = (atomid_connections[2],atomid,atomid_connections[3],atomid_connections[4]) # 2-3-4
            push!(improper_array,imp)

            imp = (atomid_connections[1],atomid,atomid_connections[2],atomid_connections[4]) # 1-2-4
            push!(improper_array,imp)

            imp = (atomid_connections[1],atomid,atomid_connections[3],atomid_connections[4]) # 1-3-4
            push!(improper_array,imp)

        end
    end
    println("Impropers determined!")
    return improper_array
end