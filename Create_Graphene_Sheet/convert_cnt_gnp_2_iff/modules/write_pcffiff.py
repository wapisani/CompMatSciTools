

def string_coeffs(coeff_list):

    string = ''
    for k in coeff_list:
        string += str(k) + ' '
    return string.strip()


def moleculefile(outfile, m):

    if m.ndihedrals != 0:
        hasdihedrals = True
    else:
        hasdihedrals = False

    try:
        if m.nimpropers != 0:
            hasimpropers = True
        else:
            hasimpropers = False
    except BaseException:
        hasimpropers = False

    f = open(outfile, 'w')

    f.write('HEADER\n\n')

    f.write('{0} atoms\n{1} bonds\n{2} angles\n{3} dihedrals\n'
            .format(m.natoms, m.nbonds, m.nangles, m.ndihedrals))

    if hasimpropers:
        f.write('{0} impropers\n'.format(m.nimpropers))

    f.write('\n')

    f.write('{0} atom types\n{1} bond types\n{2} angle types\n{3} dihedral types\n'
            .format(m.natomtypes, m.nbondtypes, m.nangletypes, m.ndihedraltypes))

    if hasimpropers:
        f.write('{0} improper types\n'.format(m.nimpropertypes))

    f.write('\n')

    if m.extra_lines is not None and m.extra_lines != '':
        f.write(m.extra_lines + '\n')

    f.write('{0}\n{1}\n{2}\n\n'.format(m.xbox_line, m.ybox_line, m.zbox_line))

    f.write('Masses\n\n')
    for i in m.masses:
        f.write('{0} {1}\n'.format(i, m.masses[i]))

    f.write('\nPair Coeffs\n\n')
    for i in m.pair_coeffs:
        string = string_coeffs(m.pair_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nBond Coeffs\n\n')
    for i in m.bond_coeffs:
        string = string_coeffs(m.bond_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nAngle Coeffs\n\n')
    for i in m.angle_coeffs:
        string = string_coeffs(m.angle_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nDihedral Coeffs\n\n')
    for i in m.dihedral_coeffs:
        string = string_coeffs(m.dihedral_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    if hasimpropers:
        f.write('\nImproper Coeffs\n\n')
        for i in m.improper_coeffs:
            string = string_coeffs(m.improper_coeffs[i])
            f.write('{0} {1}\n'.format(i, string))

    f.write('\nBondBond Coeffs\n\n')
    for i in m.bondbond_coeffs:
        string = string_coeffs(m.bondbond_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nBondAngle Coeffs\n\n')
    for i in m.bondangle_coeffs:
        string = string_coeffs(m.bondangle_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    if hasimpropers:
        f.write('\nAngleAngle Coeffs\n\n')
        for i in m.angleangle_coeffs:
            string = string_coeffs(m.angleangle_coeffs[i])
            f.write('{0} {1}\n'.format(i, string))

    f.write('\nAngleAngleTorsion Coeffs\n\n')
    for i in m.angleangletors_coeffs:
        string = string_coeffs(m.angleangletors_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nEndBondTorsion Coeffs\n\n')
    for i in m.endbondtorsion_coeffs:
        string = string_coeffs(m.endbondtorsion_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nMiddleBondTorsion Coeffs\n\n')
    for i in m.midbondtorsion_coeffs:
        string = string_coeffs(m.midbondtorsion_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nBondBond13 Coeffs\n\n')
    for i in m.bondbond13_coeffs:
        string = string_coeffs(m.bondbond13_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nAngleTorsion Coeffs\n\n')
    for i in m.angletorsion_coeffs:
        string = string_coeffs(m.angletorsion_coeffs[i])
        f.write('{0} {1}\n'.format(i, string))

    f.write('\nAtoms\n\n')
    for i in m.atoms:
        a = m.atoms[i]
        try:
            f.write('{0} {1} {2} {3} {4} {5} {6}\n'
                    .format(i, a.typenum, a.type, a.charge, a.x, a.y, a.z))
        except BaseException:
            f.write('{0} {1} {2} {3} {4} {5} {6}\n'
                    .format(i, 1, a.type, a.charge, a.x, a.y, a.z))

    # if m.velocities != {}:
    #  f.write('\nVelocities\n\n')
    #  for i in m.velocities:
    #    v = m.velocities[i]
    #    f.write('{0} {1} {2} {3}\n'\
    #      .format(i,v[0],v[1],v[2]))

    f.write('\nBonds\n\n')
    for i in m.bonds:
        b = m.bonds[i]
        f.write('{0} {1} {2} {3}\n'
                .format(i, b.type, b.atomids[0], b.atomids[1]))

    f.write('\nAngles\n\n')
    for i in m.angles:
        c = m.angles[i]
        f.write('{0} {1} {2} {3} {4}\n'
                .format(i, c.type, c.atomids[0], c.atomids[1], c.atomids[2]))

    if hasdihedrals:
        f.write('\nDihedrals\n\n')
        for i in m.dihedrals:
            d = m.dihedrals[i]
            f.write(
                '{0} {1} {2} {3} {4} {5}\n' .format(
                    i,
                    d.type,
                    d.atomids[0],
                    d.atomids[1],
                    d.atomids[2],
                    d.atomids[3]))

    if hasimpropers:
        f.write('\nImpropers\n\n')
        for i in m.impropers:
            e = m.impropers[i]
            f.write(
                '{0} {1} {2} {3} {4} {5}\n' .format(
                    i,
                    e.type,
                    e.atomids[0],
                    e.atomids[1],
                    e.atomids[2],
                    e.atomids[3]))

    f.close()
