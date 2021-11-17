#!/usr/bin/env python

# Program written by Benjamin Jensen
# Revision 1.0
# January 7th 2013
# Michigan Technological University
# 1400 Townsend Dr.
# Houghton, MI 49913


def angles(bonds):  # returns a set of tuples
    print("Determining angles... ")
    angle_set = set([])
    count = 1
    for atom1, atom2 in bonds:
        atom12set = set([atom1, atom2])
        for atom3, atom4 in bonds:
            atom34set = set([atom3, atom4])
            atom1234set = set([atom1, atom2, atom3, atom4])
            if len(atom1234set) == 3:
                centeratom = list(atom12set - (atom12set ^ atom34set))
                outsideatoms = list(atom12set ^ atom34set)
                angle = [outsideatoms[0], centeratom[0], outsideatoms[1]]
                if tuple(angle) not in angle_set:
                    angle.reverse()
                    angle_set.add(tuple(angle))
                    print(f"{count} {angle[0]} {angle[1]} {angle[2]}")
                    count += 1
    print("Angles determined!")
    return angle_set
