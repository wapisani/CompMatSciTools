
import modules.angles as angles
import modules.dihedrals as dihedrals


class Bond:
    pass


class Angle:
    pass


class Dihedral:
    pass


def dist(t=(), s=(), l=None):
    zdist = abs(t[2] - s[2])
    if zdist > 0.5 * l:
        zdist = l - zdist
    square = (t[0] - s[0])**2 + (t[1] - s[1])**2 + zdist**2
    distance = square**0.5
    return distance


def unoverlap(pset, safety_dist, l):

    delete = set([])

    for p1 in pset:
        for p2 in pset:
            distance = dist(p1, p2, l)
            if distance <= safety_dist\
                    and p1 != p2\
                    and p1 not in delete\
                    and p2 not in delete:
                delete.add(p1)

    new_set = pset.difference(delete)
    return new_set


def setify(items):  # takes a list of bonds/angles and makes them a "set" of "tuples"
    my_set = set([])
    for item in items.values():
        my_set.add(tuple(item.atomids))
    return my_set


def findbonds(a, maxdist, l):

    bonded = set([])
    bond_dict = {}

    b_id = 1
    for a1 in a:
        for a2 in a:
            if a1 != a2:
                p1 = (a[a1].x, a[a1].y, a[a1].z)
                p2 = (a[a2].x, a[a2].y, a[a2].z)
                distance = dist(p1, p2, l)
                if distance <= maxdist\
                        and (a1, a2) not in bonded\
                        and (a2, a1) not in bonded:

                    bonded.add((a1, a2))

                    b = Bond()
                    b.type = 1
                    b.atomids = [a1, a2]
                    bond_dict[b_id] = b
                    b_id += 1

    return bond_dict


def findangles(mbonds):

    c = {}

    # finds angles based on connectivity
    bond_set = setify(mbonds)
    angle_set = angles.angles(bond_set)

    # id counter
    group_id = 1

    # organize info into a dictionary of Angle objects
    for group in angle_set:
        an = Angle()
        an.type = 1  # dummy type
        an.atom1id = group[0]
        an.atom2id = group[1]
        an.atom3id = group[2]
        an.atomids = [an.atom1id, an.atom2id, an.atom3id]
        c[group_id] = an
        group_id = group_id + 1

    return c


def finddihedrals(mbonds, mangles):

    d = {}

    # finds angles based on connectivity
    bond_set = setify(mbonds)
    angle_set = setify(mangles)

    # finds dihedrals based on connectivity and previously found angles
    dihedral_set = dihedrals.dihedrals(bond_set, angle_set)

    group_id = 1

    # organize info into a dictionary of Dihedral objects
    for group in dihedral_set:
        di = Dihedral()
        di.type = 1  # dummy type
        di.atom1id = group[0]
        di.atom2id = group[1]
        di.atom3id = group[2]
        di.atom4id = group[3]
        di.atomids = [di.atom1id, di.atom2id, di.atom3id, di.atom4id]
        d[group_id] = di
        group_id = group_id + 1

    return d
