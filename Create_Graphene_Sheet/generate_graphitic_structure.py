#This python script will creat a single layer of graphite (graphene).
#Program written by Choi Sorayot Chinkanjanarot
#Created on Aug 18,2015
#Add bonds creating Aug 20,2015
#Add bonds creating over the boundary and also the angles Aug 21, 2015
#Add dihedral creating Aug 24,2015
#Add coefficient Aug 24,2015

# Re-written by Will Pisani on April 10, 2020


print('\n')
print('Graphene creator')

import sys, os
import angles
import dihedrals
import read
import write
import numpy as np
from numba import jit
import time

program_start_time = time.perf_counter()

os.chdir(r"your/directory/here")

# I'm not sure if this jit really helps or not
@jit
def dist3d(x1, y1, z1, x2, y2, z2):
     dist = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
     return dist

class Bond: pass #Class for creating bond's objects
class Angle: pass #Class for creating angle's objects
class Dihedral: pass #Class for creating dihedral's object

# Set input filename
input_filename = r"GNP_lattice_x20.01_y35.0175_test.dat"
# input_filename = r"GNP_lattice_x2.01_y3.498405_test.dat"

# Read in data structure from data file
m = read.Molecule_File(input_filename)

dxx = m.xbox_line.split() 
dyy = m.ybox_line.split() 
dxxx = 2*float(dxx[1])
dyyy = 2*float(dyy[1])
#print dxxx, dyyy
dx = [dxxx, 0, -dxxx]
dy = [dyyy, 0, -dyyy]

#Add coefficients
ncoef = 1
pc = m.pair_coeffs
bc = m.bond_coeffs
nc = m.angle_coeffs
dc = m.dihedral_coeffs
pc[ncoef] = [0.0700, 3.5500]
bc[ncoef] = [469.00, 1.400]
nc[ncoef] = [63.00, 120.00]
dc[ncoef] = [0, 7.25, 0, 0]
m.natomtypes = 1
m.nbondtypes = 1
m.nangletypes = 1
m.ndihedraltypes = 1

def update_numbers(molecule_structure):
    """This function will update the number of atoms, bonds, angle, dihedrals, impropers,
    bond types, angle types, dihedral types, improper types, etc. of a Molecule_File
    data structure."""
    mol = molecule_structure
    mol.natoms = len(mol.atoms)
    mol.nbonds = len(mol.bonds)
    mol.nangles = len(mol.angles)
    mol.ndihedrals = len(mol.dihedrals)
    mol.nbondtypes = len(mol.bond_coeffs)
    mol.nangletypes = len(mol.angle_coeffs)
    mol.ndihedraltypes = len(mol.dihedral_coeffs)
    
    return mol

# Create a new function to unpack the atom ids and coordinates into a 
# 2D numpy array
def unpack_atom(molecule_structure):
    atoms = molecule_structure.atoms
    atomid = []
    atomx = []
    atomy = []
    atomz = []
    for atom in atoms:
        atomid.append(atom)
        atomx.append(atoms[atom].x)
        atomy.append(atoms[atom].y)
        atomz.append(atoms[atom].z)

    return atomid,atomx,atomy,atomz

# These functions (convert_to_numpy and compute_distances) are thus far 
# unused and incomplete. They're meant to replace createBonds
def convert_to_numpy(atomid,atomx,atomy,atomz):
    atom_pos = np.zeros((len(atomx),3))
    for index,value in enumerate(atomx):
        atom_pos[index,0] = atomx[index]
        atom_pos[index,1] = atomy[index]
        atom_pos[index,2] = atomz[index]
    return atom_pos

def compute_distances(atom_pos,dx,dy):
    non_periodic_dist = np.sqrt(np.sum((atom_pos[0:-2] - atom_pos[1:-1])**2,axis=1))

# Define new bond creation function using the atominfo array
def createBonds(atomid,atomx,atomy,atomz,dx,dy):
    bond_list = [] # Store each bond as a 2-element list within the list
    bond_count = 1 # Counter used for print information
    for index1,atom1 in enumerate(atomid):
        x1 = atomx[index1]
        y1 = atomy[index1]
        z1 = atomz[index1]
        for index2,atom2 in enumerate(atomid):
            if atom1 != atom2:
                for delx in dx:
                    for dely in dy:
                        x2 = atomx[index2] + delx
                        y2 = atomy[index2] + dely
                        z2 = atomz[index2]
                        distance = dist3d(x1,y1,z1,x2,y2,z2)
                        if distance < 1.425:
                            if [atom1,atom2] not in bond_list or [atom2,atom1] not in bond_list:
                                
                                bond_list.append([atom1,atom2])
                                print(f"{bond_count} {atom1} {atom2}")
                                bond_count += 1
                            
    return bond_list


def GetAngles(bond_list):
    angle_set = angles.angles(bond_list)
    angle_list = [list(angle) for angle in angle_set]

    return angle_list

def GetDihedrals(bond_list,angle_list):
    dihedral_set = dihedrals.dihedrals(bond_list,angle_list)
    dihedral_list = [list(dihedral) for dihedral in dihedral_set]

    return dihedral_list

bond_start_time = time.perf_counter()

atomid,atomx,atomy,atomz = unpack_atom(m)
bond_list = createBonds(atomid,atomx,atomy,atomz,dx,dy)



# Add bonds to data structure
for index,bond in enumerate(bond_list):
    index += 1 # Index in python starts at 0
    atom1 = bond[0]
    atom2 = bond[1]
    b = Bond()
    b.atomids = [atom1,atom2]
    b.type = 1
    m.bonds[index] = b

bond_stop_time = time.perf_counter() - bond_start_time
print(f"Time taken to get a list of all bonds: {bond_stop_time} seconds\n")

angle_start_time = time.perf_counter()

angle_list = GetAngles(bond_list)



# Add angles to data structure
for index,angle in enumerate(angle_list):
    index += 1 # Index in python starts at 0
    atom1 = angle[0]
    atom2 = angle[1]
    atom3 = angle[2]
    a = Angle()
    a.atomids = [atom1,atom2,atom3]
    a.type = 1
    m.angles[index] = a
    
angle_stop_time = time.perf_counter() - angle_start_time
print(f"Time taken to get a list of angles: {angle_stop_time} seconds\n")
  
dihedral_start_time = time.perf_counter()

dihedral_list = GetDihedrals(bond_list,angle_list)



# Add dihedrals to data structure
for index,dihedral in enumerate(dihedral_list):
    index += 1 # Index in python starts at 0
    atom1 = dihedral[0]
    atom2 = dihedral[1]
    atom3 = dihedral[2]
    atom4 = dihedral[3]
    d = Dihedral()
    d.atomids = [atom1,atom2,atom3,atom4]
    d.type = 1
    m.dihedrals[index] = d

dihedral_stop_time = time.perf_counter() - dihedral_start_time
print(f"Time taken to get a list of dihedrals: {dihedral_stop_time} seconds\n")

# Update all nbonds, nangles, etc. in data structure    
m = update_numbers(m)

noa = m.natoms
nob = m.nbonds
non = m.nangles
nod = m.ndihedrals

# Set output filename
output_filename = input_filename[:-4] + "_opls.dat"

# Create a new data file
write.moleculefile(output_filename,m)


print("____________________________________________________________________\n")
print(f"Data file written to {output_filename} in directory {os.getcwd()}\n")
print("#atoms : ", noa)
print("#bonds : ", nob)
print("#angles : ", non)
print("#dihedrals : ", nod)

program_time_taken = time.perf_counter() - program_start_time
print(f"\nTotal time taken: {program_time_taken} seconds")
print("Timing breakdown for each section")
print(f"Bond creation: {bond_stop_time} seconds")
print(f"Angle creation: {angle_stop_time} seconds")
print(f"Dihedral creation: {dihedral_stop_time} seconds")

