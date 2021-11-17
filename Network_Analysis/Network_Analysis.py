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
Created on Mon Jul  6 17:21:43 2020

@author: William A. Pisani, Ph.D.

This script will read in a LAMMPS PCFF-IFF data file and load it as a network 
structure for analysis. Please be warned that this can take a while for a large system
(tens of thousands of atoms). Be sure to read through this entire script before running.

Future improvements:
Implement as a Jupyter notebook
"""

# Import necessary libraries
import os
import read_PCFF
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
import matplotlib.animation as animation

# Define functions here
def getChainCoords(chainNodes,boxLengths):
    """

    Parameters
    ----------
    chainNodes : list
        List of atom ids making up a chain of atoms.
    boxLengths : tuple
        Tuple of the simulation box length in order of (xlen, ylen, zlen)
    Returns
    -------
    uwrap_xpos : list
        A list of unwrapped x-positions in order from chain end to chain end.
    uwrap_ypos : list
        A list of unwrapped y-positions in order from chain end to chain end.
    uwrap_zpos : list
        A list of unwrapped z-positions in order from chain end to chain end.
    wrap_xpos : list
        A list of wrapped x-positions in order from chain end to chain end.
    wrap_ypos : list
        A list of wrapped y-positions in order from chain end to chain end.
    wrap_zpos : list
        A list of wrapped z-positions in order from chain end to chain end.

    """
    xlen = boxLengths[0]
    ylen = boxLengths[1]
    zlen = boxLengths[2]
    
    uwrap_xpos = []
    uwrap_ypos = []
    uwrap_zpos = []
    wrap_xpos = []
    wrap_ypos = []
    wrap_zpos = []
    
    for atom in chainNodes:
        x = m.atoms[atom].x
        y = m.atoms[atom].y
        z = m.atoms[atom].z
        wrap_xpos.append(x)
        wrap_ypos.append(y)
        wrap_zpos.append(z)
        
        # Unwrap coordinates and save them in array
        ux = x + m.atoms[atom].ix * xlen
        uy = y + m.atoms[atom].iy * ylen
        uz = z + m.atoms[atom].iz * zlen
        
        uwrap_xpos.append(ux)
        uwrap_ypos.append(uy)
        uwrap_zpos.append(uz)
        
    return uwrap_xpos, uwrap_ypos, uwrap_zpos, wrap_xpos, wrap_ypos, wrap_zpos

# Change to directory
directory = r"Your_path_here"
os.chdir(directory)
data_file = r"PA6IFFPlyEq.dat"

print(f"Analysis of the polymer network of {data_file}...")
print("Reading in the data file now ...")
# Load data file into data structure
m = read_PCFF.Molecule_File(data_file)
print("Done")
print("Now creating directed graph of the specified bonds.")

# Initialize directed graph
G=nx.DiGraph()

# Iterate over bonds and add them to the directed graph
# Specify the bond types you do not want to add to the graph.
# I recommend specifying any bonds that are not part of the molecule backbone
# like hydrogen atoms (unless that hydrogen atom is at the end of the molecule backbone
# like in the case of an alcohol).
# N6GNP bond types not used [4,5,8,11,12,14,16,17,18]
####### You will need to change these types to match your own system!!! ###############
types_not_used = [4,5,8,11,12,14,16,17,18]
for bond in m.bonds:
    atom1 = m.bonds[bond].atomids[0]
    atom2 = m.bonds[bond].atomids[1]
    bond_type = m.bonds[bond].type
    if bond_type not in types_not_used:
        G.add_edge(atom1,atom2,type=bond_type)
    
print("Finished creating directed graph.\nNow traversing the graph to find the chains...")
# Traverse the network to get a list of chains
# Original code by Arkeen at https://stackoverflow.com/questions/55711945/networkx-getting-all-possible-paths-in-dag
# From the networkx documentation, "The node in-degree is the number of edges pointing in to the node."
# and "The node out_degree is the number of edges pointing out of the node."
roots = []
leaves = []
all_paths = []
for node in G.nodes:
  if G.in_degree(node) == 0: # it's a root
    roots.append(node)
  elif G.out_degree(node) == 0: # it's a leaf
    leaves.append(node)

for root in roots:
    for leaf in leaves:
        for path in nx.all_simple_paths(G, root, leaf):
            all_paths.append(path)            
            # print(path)

# Now get all simple cycles from the graph and add to paths
# Unfortunately, the above code doesn't work for loops/cycles, but the below code
# will find them and add them to the all_paths list. The existence of cycles is
# why the code didn't match the number of clusters found by Ovito.
# Now it should match the number given by Ovito after adding the number of 
# reinforcement molecules.
print("Done!")
print("Now finding all loops...")
loop_counter = 0
loop_paths = []
for cycle in nx.simple_cycles(G):
    loop_paths.append(cycle)
    all_paths.append(cycle)
    loop_counter += 1
    
print("Done!")

# Now get MW of the cycles
# Iterate over the chains and get the molecular weight of each chain
cycles_MW_array = []
for chain in loop_paths:
    mol_weight = 0
    for atom in chain:
        type = m.atoms[atom].type
        mol_weight += m.masses[type]
    cycles_MW_array.append(mol_weight)
# Print the number of paths (chains) and the polymerization percentage based on this

print(f"\nNumber of rings/loops/cycles: {loop_counter}")
print(f"Avg MW of rings: {np.average(cycles_MW_array):.2f} +- {np.std(cycles_MW_array):.4f} g/mol")
# Get the number of nitrogen atoms, this will give the number of monomers the sims started with
nitrogen_count = 0
for atom in m.atoms:
    atomtype = m.atoms[atom].type
    if atomtype == 7 or atomtype == 12:
        nitrogen_count += 1
poly_percent = (nitrogen_count - len(all_paths))/nitrogen_count*100
print(f"\nNumber of chains: {len(all_paths)}")
print(f"Extent of polymerization: {poly_percent:.2f}%")
# nx.draw_networkx(G, pos=nx.circular_layout(G))
# plt.show()

# Get average length of chains
chain_length_array = []
for chain in all_paths:
    chain_length_array.append(len(chain))
    
avg_chain_length = np.average(chain_length_array)
shortest_chain_length = min(chain_length_array)
shortest_chain_index = chain_length_array.index(min(chain_length_array))
shortest_chain_nodes = all_paths[shortest_chain_index]
longest_chain_length = max(chain_length_array)
longest_chain_index = chain_length_array.index(max(chain_length_array))
longest_chain_nodes = all_paths[longest_chain_index]

print("\nThe hydrogen and oxygen (in the c-o double bond) are not included in the following atom counts.")
print("This is due to the chain counting mechanism. Including any atoms not in the backbone will result in extra chains that don't actually exist.")

print(f"\nAverage chain length: {avg_chain_length:.2f} atoms")
print(f"Shortest chain length: {shortest_chain_length} atoms")
print(f"Longest chain length: {longest_chain_length} atoms")

# Figure out how many molecules are in each chain, number of molecules will change 
# based on the ends
num_mol_array = []
for chain in all_paths:
    mol_count = 0
    for atom in chain:
        a = m.atoms[atom]
        if a.type == 7 or a.type == 12:
            mol_count += 1
    num_mol_array.append(mol_count)
    
print(f"\nAverage chain length: {np.average(num_mol_array):.1f} molecules")
print(f"Shortest chain length: {min(num_mol_array):.1f} molecules")
print(f"Longest chain length: {max(num_mol_array):.1f} molecules")

# Make a histogram of the chain lengths
# print("Histogram needs some work")
# hist,bin_edges = np.histogram(chain_length_array,bins=10)
# fig_hist, ax_hist = plt.subplots()
# ax_hist.bar(bin_edges[:-1],hist,width = 5, color='#0504aa',alpha=0.7)
# ax_hist.axis(xmin=0,xmax=max(bin_edges))
# ax_hist.grid(axis='y', alpha=0.75)
# ax_hist.set_xlabel('Chain length, atoms',fontsize=15)
# ax_hist.set_ylabel('Frequency',fontsize=15)
# # fig_hist.set_xticks(fontsize=15)
# # fig_hist.set_yticks(fontsize=15)
# ax_hist.set_ylabel('Number of chains',fontsize=15)
# rect_heights = ax_hist.patches

# for rect, label in zip(rect_heights, bin_edges[:-1]):
#     height = rect.get_height()
#     ax_hist.text(rect.get_x() + rect.get_width() / 2, height + 1, int(label),
#             ha='center', va='bottom')


# Now iterate over the chains and get the physical distance between connected atoms
# in each chain. This will give the physical length of the chains.
chain_distance_array = []
xlen = m.xhi - m.xlo
ylen = m.yhi - m.ylo
zlen = m.zhi - m.zlo

boxLengths = (xlen,ylen,zlen)

# The delx, dely, delz are to figure out distances of bonds straddling the periodic boundary
delx = [-xlen/2,0,xlen/2]
dely = [-ylen/2,0,ylen/2]
delz = [-zlen/2,0,zlen/2]

for chain in all_paths:
    distance = 0
    for index,atom in enumerate(chain[:-1]):
        atom1 = m.atoms[atom]
        atom2 = m.atoms[chain[index+1]]
        x1 = atom1.x
        y1 = atom1.y
        z1 = atom1.z
        
        break_flag = False
        for dx in delx:
            for dy in dely:
                for dz in delz:
                   x2 = atom2.x + dx
                   y2 = atom2.y + dy
                   z2 = atom2.z + dz
                   current_pair_dist = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
                   # Only add computed distance between atoms if it's less than a specified number
                   # This prevent extra distance being added when one atom is on one side
                   # of the box and the second atom is on the other side of the box.
                   if current_pair_dist < 10:
                       distance += current_pair_dist
                       break_flag = True
                       break
                if break_flag:
                    break
            if break_flag:
                break
    chain_distance_array.append(distance)

print("\nThe distance covered by the chains was also computed. The distance between each connected atom was calculated and summed for each chain.")
print(f"\nAverage chain length: {np.average(chain_distance_array):.2f} angstroms")
print(f"Shortest chain length: {min(chain_distance_array):.2f} angstroms")
print(f"Longest chain length: {max(chain_distance_array):.2f} angstroms")

# Iterate over chains to get vectors between each bonded pair of atoms
# This will hopefully be used for determining inflection points
chain_vector_array = []
for chain in all_paths:
    individual_chain_vector_array = []
    for index,atom in enumerate(chain[:-1]):
        atom1 = m.atoms[atom]
        atom2 = m.atoms[chain[index+1]]
        x1 = atom1.x
        y1 = atom1.y
        z1 = atom1.z
        x2 = atom2.x
        y2 = atom2.y
        z2 = atom2.z
        individual_chain_vector_array.append((x2 - x1, y2 - y1, z2 - z1))
    chain_vector_array.append(individual_chain_vector_array)
        
# Iterate over the chains and get the molecular weight of each chain
mol_weight_array = []
for chain in all_paths:
    mol_weight = 0
    for atom in chain:
        type = m.atoms[atom].type
        mol_weight += m.masses[type]
    mol_weight_array.append(mol_weight)

# Get number average molecular weight
# https://chem.libretexts.org/Bookshelves/Analytical_Chemistry/Book%3A_Physical_Methods_in_Chemistry_and_Nano_Science_(Barron)/02%3A_Physical_and_Thermal_Analysis/2.02%3A_Molecular_Weight_Determination#:~:text=Molecular%20Weight%20Calculations,-Number%20average%20of&text=Number%20average%20of%20molecular%20weight%20is%20the%20total%20weight%20of,n)%20is%20given%20by%202.2.
total_mass = 0
weight_avg_numerator = 0
for mass, num in Counter(mol_weight_array).most_common():
    total_mass += mass * num
    weight_avg_numerator += num * mass ** 2

num_avg_mw = total_mass/len(mol_weight_array)
weight_avg_mw = weight_avg_numerator/total_mass

print(f"\nNumber average molecular weight: {num_avg_mw:.2f} g/mol")    
print(f"Weight average molecular weight: {weight_avg_mw:.2f} g/mol")

# Sort and then print the molecular weights
print("\nTop 10 Chains")
print("Chain # Molecular Weight (g/mol)")
print("---------------------------------")
count = 1
sorted_mw_array = sorted(mol_weight_array,reverse=True)
for mw in sorted_mw_array:
    if count > 10:
        break
    else:
        print(f" {count}\t\t{mw:.2f}")
        count += 1

print("\n\nMolecular Weight (g/mol)       Number of chains")
print("-----------------------------------------------------")
counted_frag = Counter(mol_weight_array)
counted_frag = {k: counted_frag[k] for k in sorted(counted_frag,reverse=True)} # Create an ordered dict
for key in counted_frag:
    print(f"\t{key:.2f}\t\t\t\t{counted_frag[key]}")

# Iterate over the chains and the atoms within those chains and store the 
# x-, y-, and z-coordinates as a list of lists
all_paths_xcoord_array = []
all_paths_ycoord_array = []
all_paths_zcoord_array = []
for chain in all_paths:
    chain_xcoord_array = []
    chain_ycoord_array = []
    chain_zcoord_array = []
    for atom in chain:
        a = m.atoms[atom]
        chain_xcoord_array.append(a.x)
        chain_ycoord_array.append(a.y)
        chain_zcoord_array.append(a.z)
    all_paths_xcoord_array.append(chain_xcoord_array)
    all_paths_ycoord_array.append(chain_ycoord_array)
    all_paths_zcoord_array.append(chain_zcoord_array)

# Iterate over these chain coordinates and determine the xmin, xmax, ymin, ymax,
# zmin, and zmax for each chain. Also determine what percentage of each dimension's box
# length the chain spans the box in 
# each dimension and store the answer as the last three positions in the following
# list of tuples:  [(xmin,xmax,ymin,ymax,zmin,zmax,xspan,yspan,zspan)]
# This is a measure that I came up with that may or may not be sound.
# It involves determining the maximum and minimum coordinates of each chain's atoms.
# So the max and min of each chain in each dimension likely comes from six different
# atoms in each chain. I'm trying to get a sense of how much each chain spans
# the box dimensions.
all_chains_coord_boundaries = []
for index,chain_coords in enumerate(all_paths_xcoord_array):
    xmin = min(all_paths_xcoord_array[index])
    xmax = max(all_paths_xcoord_array[index])
    ymin = min(all_paths_ycoord_array[index])
    ymax = max(all_paths_ycoord_array[index])
    zmin = min(all_paths_zcoord_array[index])
    zmax = max(all_paths_zcoord_array[index])
    chain_xlen = xmax - xmin
    chain_ylen = ymax - ymin
    chain_zlen = zmax - zmin
    xspan = chain_xlen/xlen
    yspan = chain_ylen/ylen
    zspan = chain_zlen/zlen
    all_chains_coord_boundaries.append((xmin,xmax,ymin,ymax,zmin,zmax,xspan,yspan,zspan))


# Plot molecular weight vs span
chains_xspan = [boundary[6] for boundary in all_chains_coord_boundaries]
chains_yspan = [boundary[7] for boundary in all_chains_coord_boundaries]
chains_zspan = [boundary[8] for boundary in all_chains_coord_boundaries]
chains_avgspan = [(boundary[6]+boundary[7]+boundary[8])/3 for boundary in all_chains_coord_boundaries]

# Get the 

fig_mwspan, ax_mwspan = plt.subplots()
# ax_mwspan.scatter(mol_weight_array,chains_xspan,marker=".",facecolor="none",color="orange",label="X")
# ax_mwspan.scatter(mol_weight_array,chains_yspan,marker="s",facecolor="none",color="blue",label="Y")
# ax_mwspan.scatter(mol_weight_array,chains_zspan,marker="p",facecolor="none",color="green",label="Z")
ax_mwspan.scatter(mol_weight_array,chains_avgspan,marker="*",color="black",s=80,label="AVG")
ax_mwspan.set_xlabel("Molecular weight, g/mol")
ax_mwspan.set_ylabel("Span")
ax_mwspan.legend()

# Get chain coordinates of the shortest chain (though there may be multiple shortest chains)
shortest_xpos, shortest_ypos, shortest_zpos,\
    shortest_wrap_xpos, shortest_wrap_ypos, shortest_wrap_zpos = getChainCoords(shortest_chain_nodes,boxLengths)


# Get chain coordinates of the longest chain (there might be multiple longest chains)
longest_xpos, longest_ypos, longest_zpos,\
    longest_wrap_xpos, longest_wrap_ypos, longest_wrap_zpos = getChainCoords(longest_chain_nodes,boxLengths)  



# Plot the shortest chain on a 3D plot, unwrapped coordinates
# fig_short = plt.figure()
# ax_short = fig_short.add_subplot(111,projection='3d')
# ax_short.plot(shortest_xpos,shortest_ypos,shortest_zpos, marker='x')#,facecolor=(0,0,0,0),edgecolors='blue')
# ax_short.set_title(f"Shortest chain is {shortest_chain_length/9:.2f} molecules")
# ax_short.set_xlabel('Unwrapped Coords, A')
# ax_short.set_ylabel('Unwrapped Coords, A')
# ax_short.set_zlabel('Unwrapped Coords, A')

# Plot the longest chain on a 3D plot, unwrapped coordinates
# fig_long = plt.figure()
# ax_long = fig_long.add_subplot(111,projection='3d')
# ax_long.plot(longest_xpos,longest_ypos,longest_zpos, marker='x')#,facecolor=(0,0,0,0),edgecolors='blue')
# ax_long.set_title(f"Longest chain is {longest_chain_length/9:.2f} molecules")
# ax_long.set_xlabel('Unwrapped Coords, A')
# ax_long.set_ylabel('Unwrapped Coords, A')
# ax_long.set_zlabel('Unwrapped Coords, A')



# Plot the shortest chain on a 3D plot, wrapped coordinates
# fig_short = plt.figure()
# ax_short = fig_short.add_subplot(111,projection='3d')
# ax_short.plot(shortest_wrap_xpos,shortest_wrap_ypos,shortest_wrap_zpos, marker='x')
# ax_short.set_title(f"Shortest chain is {shortest_chain_length/9:.2f} molecules")
# ax_short.set_xlabel('Wrapped Coords, A')
# ax_short.set_ylabel('Wrapped Coords, A')
# ax_short.set_zlabel('Wrapped Coords, A')

# Plot the longest chain on a 3D plot, wrapped coordinates
fig_long = plt.figure()
ax_long = fig_long.add_subplot(111,projection='3d')
ax_long.scatter(longest_wrap_xpos,longest_wrap_ypos,longest_wrap_zpos, marker='x',label="Chain",s=70)
ax_long.scatter(longest_wrap_xpos[0],longest_wrap_ypos[0],longest_wrap_zpos[0], marker='s',s=80,color="green",label="Start of chain")
ax_long.scatter(longest_wrap_xpos[-1],longest_wrap_ypos[-1],longest_wrap_zpos[-1], marker='^',s=80,color="red",label="End of chain")
ax_long.set_title(f"{max(num_mol_array):.1f} molecules long, {max(chain_distance_array):.2f} angstroms long")
# ax_long.axis(xmin=m.xlo,xmax=m.xhi,ymin=m.ylo,ymax=m.yhi,zmin=m.zlo,zmax=m.zhi)
ax_long.set_xlabel('Wrapped X Coords, A')
ax_long.set_ylabel('Wrapped Y Coords, A')
ax_long.set_zlabel('Wrapped Z Coords, A')
plt.tight_layout()
ax_long.legend(loc='lower left')

# For rotating the figure, need to use %matplotlib qt in the console before running, or set backend to tkinter
for angle in range(0,360):
    ax_long.view_init(30,angle)
    plt.draw()
    plt.pause(0.001)
    
# Attempt to save a GIF of a rotation animation of the longest chain 3D plot
# def rotate_plot(self):
#     for angle in range(0,360):
#         ax_long.view_init(30,angle)
#         plt.draw()
#         plt.pause(0.001)

# rotate_anim = animation.FuncAnimation(fig_long,rotate_plot,frames=30,repeat=False)
# plt.show()
# writergif = animation.PillowWriter(fps=30)
# rotate_anim.save(data_file.split('.')[0]+'.gif',writer=writergif)
# rotate_anim.save(data_file.split('.')[0]+'.mp4',writer=FFwriter)

# Plot directed graph
# pos = nx.layout.spring_layout(G)
# M = G.number_of_edges()
# edge_colors = range(2, M + 2)
# nodes = nx.draw_networkx_nodes(G, pos, node_color='blue')
# edges = nx.draw_networkx_edges(G, pos, arrowstyle='->',
#                                arrowsize=10, edge_color=edge_colors,
#                                edge_cmap=plt.cm.Blues, width=2)
# pc = mpl.collections.PatchCollection(edges, cmap=plt.cm.Blues)
# pc.set_array(edge_colors)
# plt.colorbar(pc)

# ax = plt.gca()
# ax.set_axis_off()
# plt.show()