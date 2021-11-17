#!/usr/bin/env python

# Program written by Benjamin Jensen
# Revision 1.0
# January 7th 2013
# Michigan Technological University
# 1400 Townsend Dr.
# Houghton, MI 49913

#


class Atom:
    pass


class Bond:
    pass  # .type .atomids = [atom1id, atom2id]


class Angle:
    pass  # .type .atomids = [atom1id, atom2id, atom3id]


class Dihedral:
    pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id]


def strip_comment(line):
    # remove comments
    end = line.find('#')
    if end >= 0:
        #comment = line[end:]
        line = line[:end]
    return line


class Molecule_File:
    def __init__(self, inmolfile):
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.angles = {}  # {angle number : angle object}
        self.dihedrals = {}  # {dihedral number : dihedral object}

        self.velocities = {}  # {atom number : tuple of velocities}

        # Parameters
        self.masses = {}  # {atom type : atom mass}
        self.pair_coeffs = {}  # {atom type : list of coeffs}
        # {bond type : list of coeffs}  {1: [340,1.5], 2: [450,1.2], ...}
        self.bond_coeffs = {}
        self.angle_coeffs = {}  # {angle type : list of coeffs}
        self.dihedral_coeffs = {}  # {dihedral type : list of coeffs}

        self.total_line = ''
        self.ttype_line = ''
        self.xbox_line = ''
        self.ybox_line = ''
        self.zbox_line = ''
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

        self.parsefile(inmolfile)
        self.inmolfile = inmolfile

    def parsefile(self, inmolfile):
        f = open(inmolfile, 'r')

        massflag = False
        coeff_flag = False
        atomflag = False
        bondflag = False
        angleflag = False
        dihedralflag = False
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
                continue
            elif 'dihedrals' in line:
                self.ndihedrals = int(line.split()[0])
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
            elif 'per atom' in line:
                self.extra_lines += line + '\n'
                continue
            elif 'xlo' in line:
                self.xbox_line = line
                continue
            elif 'ylo' in line:
                self.ybox_line = line
                continue
            elif 'zlo' in line:
                self.zbox_line = line
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
            elif line == 'Atoms':
                atomflag = True

                atomstyle = str(whole_line.split('#')[-1])

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
                idcoeffs = []
                for i in line[1:]:
                    if '.' not in i:
                        idcoeffs.append(int(i))
                    else:
                        idcoeffs.append(float(i))
                        coeffs[id] = idcoeffs

            if atomflag:
                line = line.split()
                atomstyle = "molecular"
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
                    typenum = int(line[1])
                    type = int(line[2])
                    x = float(line[3])
                    y = float(line[4])
                    z = float(line[5])
                    a = Atom()
                    a.typenum = typenum
                    a.type = type
                    a.x = x
                    a.y = y
                    a.z = z
                    self.atoms[id] = a

                elif atomstyle == "full":
                    id = int(line[0])
                    typenum = int(line[1])
                    type = int(line[2])
                    charge = float(line[3])
                    x = float(line[4])
                    y = float(line[5])
                    z = float(line[6])
                    a = Atom()
                    a.typenum = typenum
                    a.type = type
                    a.charge = charge
                    a.x = x
                    a.y = y
                    a.z = z
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

            elif velocityflag:
                line = line.split()
                id = int(line[0])
                vx = float(line[1])
                vy = float(line[2])
                vz = float(line[3])
                self.velocities[id] = (vx, vy, vz)

        f.close()
