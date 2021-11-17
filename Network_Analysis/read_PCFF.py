#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Benjamin Jensen, with modifications by Will Pisani August 14, 2019
Revision 2.0  
January 30th, 2019
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49913
"""

class Atom:
    """Attributes are charge, molecule id (molid), type, x, y, and z coordinates
    and image flags ix, iy, iz
    """
    pass


class Bond:
    pass  # .type .atomids = [atom1id, atom2id]


class Angle:
    pass  # .type .atomids = [atom1id, atom2id, atom3id]


class Dihedral:
    pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id]


class Improper:
    pass  # .type .atomids = [atom1,atom2,atom3,atom4]

class Coeff_class:
    pass  # .type .coeffs = []

def update_numbers(mol):
    """This function will update the number of atoms, bonds, angle, dihedrals, impropers,
    bond types, angle types, dihedral types, improper types, etc. of a Molecule_File
    data structure."""
    mol.natoms = len(mol.atoms)
    mol.nbonds = len(mol.bonds)
    mol.nangles = len(mol.angles)
    mol.ndihedrals = len(mol.dihedrals)
    mol.nimpropers = len(mol.impropers)
    mol.nbondtypes = len(mol.bond_coeffs)
    mol.nangletypes = len(mol.angle_coeffs)
    mol.ndihedraltypes = len(mol.dihedral_coeffs)
    mol.nimpropertypes = len(mol.improper_coeffs)
    
    
    return mol


def strip_comment(line):
    # remove comments
    end = line.find('#')
    if end >= 0:
#        comment = line[end:]
        line = line[:end]
    return line


class Molecule_File:
    def __init__(self, inmolfile):
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.angles = {}  # {angle number : angle object}
        self.dihedrals = {}  # {dihedral number : dihedral object}
        self.impropers = {}  # {improper number : improper object}
        self.molecules = [] # [moleculeid1, moleculeid2, ..., moleculeidN]
        self.velocities = {}  # {atom number : tuple of velocities}

        # Parameters
        self.masses = {}  # {atom type : atom mass}
        self.pair_coeffs = {}  # {atom type : list of coeffs}
        # {bond type : list of coeffs}  {1: [340,1.5], 2: [450,1.2], ...}
        self.bond_coeffs = {}
        self.angle_coeffs = {}  # {angle type : list of coeffs}
        self.dihedral_coeffs = {}  # {dihedral type : list of coeffs}
        self.improper_coeffs = {}  # {improper type : list of coeffs}
        self.bondbond_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.bondangle_coeffs = {} # {cross-term number : list of coeffs}, angles
        self.angleangle_coeffs = {} # {cross-term number : list of coeffs}, impropers
        self.angleangletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.endbondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.middlebondtorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.bondbond13_coeffs = {} # {cross-term number : list of coeffs}, dihedrals
        self.angletorsion_coeffs = {} # {cross-term number : list of coeffs}, dihedrals

        self.total_line = ''
        self.ttype_line = ''
        self.xbox_line = ''
        self.xlo = 0
        self.xhi = 0
        self.ybox_line = ''
        self.ylo = 0
        self.yhi = 0
        self.zbox_line = ''
        self.zlo = 0
        self.zhi = 0
        self.extra_lines = ''

        self.total = 0
        self.natoms = 0
        self.natomtypes = 0
        self.nbonds = 0
        self.nbondtypes = 0
        self.nangles = 0
        self.nangletypes = 0
        self.ndihedrals = 0
        self.ndihedraltypes = 0
        self.nimpropers = 0
        self.nimpropertypes = 0
        self.nbondbond = 0
        self.nbondangle = 0
        self.nangleangle = 0
        self.nangleangletorsion = 0
        self.nendbondtorsion = 0
        self.nmiddlebondtorsion = 0
        self.nbondbond13 = 0
        self.nangletorsion = 0

        self.parsefile(inmolfile)
        self.inmolfile = inmolfile

    def parsefile(self, inmolfile):
        with open(inmolfile, 'r') as f:

            massflag = False
            coeff_flag = False
            atomflag = False
            bondflag = False
            angleflag = False
            dihedralflag = False
            improperflag = False
            velocityflag = False
    
            skip = 0
            for whole_line in f:
                # skip comment lines
                skip -= 1
                if skip >= 0:
                    continue
    
                # remove comments
                line = strip_comment(whole_line)
                line = line.strip()
    
                # begining of a section, flag the start and skip one line
                if line == '':
                    massflag = False
                    coeff_flag = False
                    atomflag = False
                    bondflag = False
                    angleflag = False
                    dihedralflag = False
                    improperflag = False
                    velocityflag = False
                elif 'atoms' in line:
                    self.total_line = line
                    self.total = int(line.split()[0])
                    self.natoms = self.total
                    continue
                elif 'bonds' in line:
                    self.nbonds = int(line.split()[0])
                    continue
                elif 'angles' in line:
                    self.nangles = int(line.split()[0])
                    self.nbondbond = int(line.split()[0])
                    self.nbondangle = int(line.split()[0])
                    continue
                elif 'dihedrals' in line:
                    self.ndihedrals = int(line.split()[0])
                    self.nangleangletorsion = int(line.split()[0])
                    self.nendbondtorsion = int(line.split()[0])
                    self.nmiddlebondtorsion = int(line.split()[0])
                    self.nbondbond13 = int(line.split()[0])
                    self.nangletorsion = int(line.split()[0])
                    continue
                elif 'impropers' in line:
                    self.nimpropers = int(line.split()[0])
                    self.nangleangle = int(line.split()[0])
                    continue
                elif 'atom types' in line:
                    self.ttype_line = line
                    self.natomtypes = int(line.split()[0])
                    continue
                elif 'bond types' in line:
                    self.nbondtypes = int(line.split()[0])
                    continue
                elif 'angle types' in line:
                    self.nangletypes = int(line.split()[0])
                    continue
                elif 'dihedral types' in line:
                    self.ndihedraltypes = int(line.split()[0])
                    continue
                elif 'improper types' in line:
                    self.nimpropertypes = int(line.split()[0])
                    continue
                elif 'per atom' in line:
                    self.extra_lines += line + '\n'
                    continue
                elif 'xlo' in line:
                    self.xbox_line = line
                    self.xlo = float(line.split()[0])
                    self.xhi = float(line.split()[1])
                    continue
                elif 'ylo' in line:
                    self.ybox_line = line
                    self.ylo = float(line.split()[0])
                    self.yhi = float(line.split()[1])
                    continue
                elif 'zlo' in line:
                    self.zbox_line = line
                    self.zlo = float(line.split()[0])
                    self.zhi = float(line.split()[1])
                    continue
                elif line == 'Masses':
                    massflag = True
                    skip = 1
                    continue
                elif line == 'Pair Coeffs':
                    coeff_flag = True
                    coeffs = self.pair_coeffs
                    skip = 1
                    continue
                elif line == 'Bond Coeffs':
                    coeff_flag = True
                    coeffs = self.bond_coeffs
                    skip = 1
                    continue
                elif line == 'Angle Coeffs':
                    coeff_flag = True
                    coeffs = self.angle_coeffs
                    skip = 1
                    continue
                elif line == 'Dihedral Coeffs':
                    coeff_flag = True
                    coeffs = self.dihedral_coeffs
                    skip = 1
                    continue
                elif line == 'Improper Coeffs':
                    coeff_flag = True
                    coeffs = self.improper_coeffs
                    skip = 1
                    continue
                elif line == 'BondBond Coeffs':
                    coeff_flag = True
                    coeffs = self.bondbond_coeffs
                    skip = 1
                    continue
                elif line == 'BondAngle Coeffs':
                    coeff_flag = True
                    coeffs = self.bondangle_coeffs
                    skip = 1
                    continue
                elif line == 'AngleAngle Coeffs':
                    coeff_flag = True
                    coeffs = self.angleangle_coeffs
                    skip = 1
                    continue
                elif line == 'AngleAngleTorsion Coeffs':
                    coeff_flag = True
                    coeffs = self.angleangletorsion_coeffs
                    skip = 1
                    continue
                elif line == 'EndBondTorsion Coeffs':
                    coeff_flag = True
                    coeffs = self.endbondtorsion_coeffs
                    skip = 1
                    continue
                elif line == 'MiddleBondTorsion Coeffs':
                    coeff_flag = True
                    coeffs = self.middlebondtorsion_coeffs
                    skip = 1
                    continue
                elif line == 'BondBond13 Coeffs':
                    coeff_flag = True
                    coeffs = self.bondbond13_coeffs
                    skip = 1
                    continue
                elif line == 'AngleTorsion Coeffs':
                    coeff_flag = True
                    coeffs = self.angletorsion_coeffs
                    skip = 1
                    continue
                elif line == 'Atoms':
                    atomflag = True
                    atomstyle = whole_line.split()[-1]
                    skip = 1
                    continue
                elif line == 'Bonds':
                    bondflag = True
                    skip = 1
                    continue
                elif line == 'Angles':
                    angleflag = True
                    skip = 1
                    continue
                elif line == 'Dihedrals':
                    dihedralflag = True
                    skip = 1
                    continue
                elif line == 'Impropers':
                    improperflag = True
                    skip = 1
                    continue
                elif line == 'Velocities':
                    velocityflag = True
                    skip = 1
                    continue
    
                if massflag:
                    line = line.split()
                    id = int(line[0])
                    mass = float(line[1])
                    self.masses[id] = mass
    
                if coeff_flag:
                    line = line.split()
                    id = int(line[0])
                    c = Coeff_class()
                    idcoeffs = []
                    bond_type = ''
                    iffr_flag = 0
                    if line[1] == 'morse' or line[1] == 'class2':
                        bond_type = line[1]
                        iffr_flag = 1
                        
                    # If the data file is in IFF-R format
                    if iffr_flag == 1:
                        
                        for i in line[2:]:
                            if '.' not in i:
                                idcoeffs.append(int(i))
                            else:
                                idcoeffs.append(float(i))
                                
                        # Update c.type with morse or class2
                        c.type = bond_type
                        
                    else: # If it is PCFF-IFF format
                        for i in line[1:]:
                            if '.' not in i and 'e' not in i:
                                idcoeffs.append(int(i))
                            else:
                                idcoeffs.append(float(i))
                        # I'm not sure what Matt/Ben envisioned for the c.type
                        # here, but their vision will only come into play for PCFF-IFF.
                        if '#' not in whole_line:
                            c.type = 'N/A'
                        else:
                            c.type = whole_line.split('#')[-1].rstrip().lstrip()
                    
                    c.coeffs = idcoeffs
                    coeffs[id] = c
    
                if atomflag:
                    line = line.split()
                    
                    if atomstyle == "charge":
                        id = int(line[0])
                        type = int(line[1])
                        charge = float(line[2])
                        x = float(line[3])
                        y = float(line[4])
                        z = float(line[5])
                        a = Atom()
                        a.type = type
                        a.charge = charge
                        a.x = x
                        a.y = y
                        a.z = z
                        self.atoms[id] = a
    
                    elif atomstyle == "molecular":
                        id = int(line[0])
                        molid = int(line[1])
                        type = int(line[2])
                        x = float(line[3])
                        y = float(line[4])
                        z = float(line[5])
                        a = Atom()
                        a.molid = molid
                        a.type = type
                        a.x = x
                        a.y = y
                        a.z = z
                        self.atoms[id] = a
    
                    elif atomstyle == "full":
                        id = int(line[0])
                        molid = int(line[1])
                        if molid not in self.molecules:
                            self.molecules.append(molid)
                        type = int(line[2])
                        charge = float(line[3])
                        x = float(line[4])
                        y = float(line[5])
                        z = float(line[6])
                        a = Atom()
                        a.molid = molid
                        a.type = type
                        a.charge = charge
                        a.x = x
                        a.y = y
                        a.z = z
                        if len(line) == 10:
                            ix = int(line[7])
                            iy = int(line[8])
                            iz = int(line[9])
                            a.ix = ix
                            a.iy = iy
                            a.iz = iz
                        self.atoms[id] = a
    
                    else:
                        print(atomstyle)
                        raise Exception('Atom Style Error')
    
                elif bondflag:
                    line = line.split()
                    id = int(line[0])
                    type = int(line[1])
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    b = Bond()
                    b.type = type
                    b.atomids = [atom1id, atom2id]
                    self.bonds[id] = b
    
                elif angleflag:
                    line = line.split()
                    id = int(line[0])
                    type = int(line[1])
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    atom3id = int(line[4])
                    c = Angle()
                    c.type = type
                    c.atomids = [atom1id, atom2id, atom3id]
                    self.angles[id] = c
    
                elif dihedralflag:
                    line = line.split()
                    id = int(line[0])
                    type = int(line[1])
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    atom3id = int(line[4])
                    atom4id = int(line[5])
                    d = Dihedral()
                    d.type = type
                    d.atomids = [atom1id, atom2id, atom3id, atom4id]
                    self.dihedrals[id] = d
    
                elif improperflag:
                    line = line.split()
                    id = int(line[0])
                    type = int(line[1])
                    atom1id = int(line[2])
                    atom2id = int(line[3])
                    atom3id = int(line[4])
                    atom4id = int(line[5])
                    i = Improper()
                    i.type = type
                    i.atomids = [atom1id, atom2id, atom3id, atom4id]
                    self.impropers[id] = i
    
                elif velocityflag:
                    line = line.split()
                    id = int(line[0])
                    vx = float(line[1])
                    vy = float(line[2])
                    vz = float(line[3])
                    self.velocities[id] = (vx, vy, vz)

        


if __name__ == '__main__':
    import os
    
#    import write_PCFF as write
#
#    # If on Windows
#    if os.name == 'nt':
#        os.chdir('..\\..\\')
#    else:  # If on Linux
#        os.chdir('../../')
#
#    file_in = 'EEC_test.dat'
#    m = Molecule_File(file_in)
#    keys_to_delete = []
#    dihedral_types_deleted = []  # The dihedral types used by these weird dihedrals
#    # aren't used by any other dihedrals and should be deleted. NOTE: The
#    # associated crossterms will also change
#
#    # Find all weird dihedrals in data file and add key to list
#    for key in m.dihedrals:
#        atom1 = m.dihedrals[key].atomids[0]
#        atom2 = m.dihedrals[key].atomids[1]
#        atom3 = m.dihedrals[key].atomids[2]
#        atom4 = m.dihedrals[key].atomids[3]
#        if atom1 == atom2 or atom1 == atom3 or atom1 == atom4 or atom2 == atom3 or atom2 == atom4 or atom3 == atom4:
#            keys_to_delete.append(key)
#
#    # Remove marked keys from dictionary
#    for key in keys_to_delete:
#        del(m.dihedrals[key])
#
#    # Renumber keys
#    d = {}
#    counter = 1
#    for key in m.dihedrals:
#        d[counter] = m.dihedrals[key]
#        counter += 1
#    m.dihedrals = d
#    # Write new data file
#    file_out = 'EEC_dihedrals.dat'
#    write.moleculefile(file_out, a)
