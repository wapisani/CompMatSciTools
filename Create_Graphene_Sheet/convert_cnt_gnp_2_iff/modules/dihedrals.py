#!/usr/bin/env python

# Program written by Benjamin Jensen
# Revision 1.0
# January 7th 2013
# Michigan Technological University
# 1400 Townsend Dr.
# Houghton, MI 49913

import modules.angles as angles


def dihedrals(bonds, nangles=False):  # returns a set of tuples
    if nangles:
        nangles = angles.angles(bonds)
    print("determining dihedrals...")
    dihedrals = set([])
    for atom1, atom2, atom3 in nangles:
        atom123set = set([atom1, atom2, atom3])

        for atom4, atom5 in bonds:
            atom45set = set([atom4, atom5])
            atom1345set = set([atom1, atom3, atom4, atom5])
            atom12345set = set([atom1, atom2, atom3, atom4, atom5])
            if len(atom12345set) == 4 and len(atom1345set) == 3:
                commonatom = list(atom45set - (atom123set ^ atom45set))
                newatom = list(atom45set - set(commonatom))
                if commonatom[0] == atom1:
                    dihedral = [newatom[0], atom1, atom2, atom3]
                else:
                    dihedral = [atom1, atom2, atom3, newatom[0]]
                if tuple(dihedral) not in dihedrals:
                    dihedral.reverse()
                    dihedrals.add(tuple(dihedral))
    print("dihedrals finished")
    return dihedrals
