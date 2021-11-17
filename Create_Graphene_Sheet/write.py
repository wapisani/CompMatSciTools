

def string_coeffs(coeff_list):

  string = ''
  for k in coeff_list:
        string += str(k) + ' '
  return string.strip()


def moleculefile(outfile,m):

  f = open(outfile,'w')

  f.write('HEADER\n\n')

  f.write('{0} atoms\n{1} bonds\n{2} angles\n{3} dihedrals\n\n{4} atom types\n{5} bond types\n{6} angle types\n{7} dihedral types\n\n'.format(m.natoms,m.nbonds,m.nangles,m.ndihedrals,m.natomtypes,m.nbondtypes,
m.nangletypes,m.ndihedraltypes))

  if m.extra_lines != None and m.extra_lines != '':
    f.write(m.extra_lines + '\n')

  f.write('{0}\n{1}\n{2}\n\n'.format(m.xbox_line,m.ybox_line,m.zbox_line))

  f.write('Masses\n\n')
  for i in m.masses:
    f.write('{0} {1}\n'.format(i,m.masses[i]))

  f.write('\nPair Coeffs\n\n')
  for i in m.pair_coeffs:
    string = string_coeffs(m.pair_coeffs[i])
    f.write('{0} {1}\n'.format(i,string))

  f.write('\nBond Coeffs\n\n')
  for i in m.bond_coeffs:
    string = string_coeffs(m.bond_coeffs[i])
    f.write('{0} {1}\n'.format(i,string))

  f.write('\nAngle Coeffs\n\n')
  for i in m.angle_coeffs:
    string = string_coeffs(m.angle_coeffs[i])
    f.write('{0} {1}\n'.format(i,string))

  f.write('\nDihedral Coeffs\n\n')
  for i in m.dihedral_coeffs:
    string = string_coeffs(m.dihedral_coeffs[i])
    f.write('{0} {1}\n'.format(i,string))

  f.write('\nAtoms\n\n')
  for i in m.atoms:
    a = m.atoms[i]
    try:
	    f.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(i,a.typenum,a.type,a.charge,a.x,a.y,a.z))
    except:
	    f.write('{0} {1} {2} {3} {4} {5}\n'\
	     .format(i,a.typenum,a.type,a.x,a.y,a.z))

  if m.velocities != {}:
    f.write('\nVelocities\n\n')
    for i in m.velocities:
      v = m.velocities[i]
      f.write('{0} {1} {2} {3}\n'\
        .format(i,v[0],v[1],v[2]))

  f.write('\nBonds\n\n')
  for i in m.bonds:
    b = m.bonds[i]
    f.write('{0} {1} {2} {3}\n'\
      .format(i,b.type,b.atomids[0],b.atomids[1]))
  
  f.write('\nAngles\n\n')
  for i in m.angles:
    c = m.angles[i]
    f.write('{0} {1} {2} {3} {4}\n'\
      .format(i,c.type,c.atomids[0],c.atomids[1],c.atomids[2]))

  f.write('\nDihedrals\n\n')
  for i in m.dihedrals:
    d = m.dihedrals[i]
    f.write('{0} {1} {2} {3} {4} {5}\n'\
      .format(i,d.type,d.atomids[0],d.atomids[1],d.atomids[2]\
      ,d.atomids[3]))

  f.close()
