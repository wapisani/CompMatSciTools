# Generating topology from a list of bonds
This tool is useful when writing fix bond/react templates and monomer data files by hand. Simply draw out your molecule, number the atoms, and write up a list of the bonds in a text file. 
The format of this text file should be:  
bond_number atom1 atom2  
1 1 2  
2 2 3  
3 2 4  
4 4 5  
etc. 

But the format of the file can also be:  
bond_number bond_type atom1 atom2  
1 1 1 2  
2 2 2 3  
3 2 3 4  
etc. 

Please see PA6IFFPly_Example_Bonds.txt as an example.

This tool requires Julia 1.5.2 or greater to run and the StatsBase package. For installation of Julia, visit the [Julia website](https://julialang.org/). The Julia files angles.jl, dihedrals.jl, impropers.jl, and write_topology.jl must be in the same directory as Get_Topology_From_Bonds.jl.

To run this tool, in the terminal or command prompt, run 
`julia Get_Topology_From_Bonds.jl your_bond_list_file_here.txt`. 
As an example, try running
`julia Get_Topology_From_Bonds.jl PA6IFFPly_Example_Bonds.txt`