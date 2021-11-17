

import sys
import modules.read_molecular as read_molecular 
import modules.write_pcffiff as write_pcffiff
import modules.locality as locality
import md_math

class Molecule: pass
class Atom:     pass
class Bond:     pass
class Angle:    pass
class Dihedral: pass

file_in = sys.argv[1]
file_out = file_in + ".pcff-iff"

m = read_molecular.Molecule_File(file_in)
a = m.atoms
b = m.bonds

cntatoms = m.natoms
cnttypes = m.natomtypes
vaid = cntatoms + 1

# Get box dims
xline = m.xbox_line.split()
yline = m.ybox_line.split()
zline = m.zbox_line.split()
lx = float(xline[1])-float(xline[0])
ly = float(yline[1])-float(yline[0])
lz = float(zline[1])-float(zline[0])

##############################
### Overwrite coefficients ###
##############################

# Masses & Pair Coeffs
for atype in range(1,cnttypes+1):
  m.masses[atype]          = 10.011150
  m.masses[atype+cnttypes] = 1.000000

  m.pair_coeffs[atype]          = ['0.06200','3.93200']
  m.pair_coeffs[atype+cnttypes] = ['0.00001','0.00001']

# Bonds
# 1: cg1-cg1, 2: cg1-cge
m.bond_coeffs	= {1:['1.42','480.0','0.0','0.0'],2:['0.65','250.0','0.0','0.0']}

# Angles
# 1: cg1-cg1-cg1, 2: cg1-cg1-cge, 3: cge-cg1-cge
m.angle_coeffs	= {1:['120.0','90.0','0.0','0.0'],2:['90.0','50.0','0.0','0.0'],3:['180.0','50.0','0.0','0.0']}

# Dihedrals
m.dihedral_coeffs = {1:['0.0','0.0','0.0','0.0','0.0','0.0']}

# Impropers
m.improper_coeffs = {1:['0.0','0.0']}

# Cross Terms
m.bondbond_coeffs       = {1:['0.0','0.0','0.0'],2:['0.0','0.0','0.0'],3:['0.0','0.0','0.0']}
m.bondangle_coeffs      = {1:['0.0','0.0','0.0','0.0'],2:['0.0','0.0','0.0','0.0'],3:['0.0','0.0','0.0','0.0']}
m.angleangle_coeffs     = {1:['0.0','0.0','0.0','0.0','0.0','0.0']}
m.angleangletors_coeffs = {1:['0.0','0.0','0.0']}
m.endbondtorsion_coeffs = {1:['0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0']}
m.midbondtorsion_coeffs = {1:['0.0','0.0','0.0','0.0']}
m.bondbond13_coeffs     = {1:['0.0','0.0','0.0']}
m.angletorsion_coeffs   = {1:['0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0']}


###############################
####### Find C-C Bonds ########
###############################

ccbonds = locality.findbonds(a,2,lz)
b.update(ccbonds)
cntbonds = len(m.bonds)
vbid = cntbonds + 1

# Find which atoms are bonded to a given atom
def findbonded(a):
  bonded = []
  for id in b:
    if a == b[id].atomids[0]:
      bonded.append(b[id].atomids[1])
    elif a == b[id].atomids[1]:
      bonded.append(b[id].atomids[0])
  return bonded


################################
### Create New Atoms & Bonds ###
################################

# Reset C-C bond types
for id in b:
  b[id].type = 1

for id in xrange(1,cntatoms+1):

  # Add charge to Carbon atom
  # Account for internal polarity
  a[id].charge = 0.200

  # Approximately place virtual atoms
  # Find plane created by nearest neighbors
  bonded = findbonded(id)
  points = [[a[id].x,a[id].y,a[id].z]]
  for id2 in bonded:
    diff = a[id].x - a[id2].x
    if diff > 0.5*lx:
	px = a[id2].x + lx
    elif diff < -0.5*lx:
	px = a[id2].x - lx
    else:
	px = a[id2].x	
    
    diff = a[id].y - a[id2].y
    if diff > 0.5*ly:
	py = a[id2].y + ly
    elif diff < -0.5*ly:
	py = a[id2].y - ly
    else:
	py = a[id2].y	
    
    diff = a[id].z - a[id2].z
    if diff > 0.5*lz:
	pz = a[id2].z + lz
    elif diff < -0.5*lz:
	pz = a[id2].z - lz
    else:
	pz = a[id2].z
	
    points.append([px,py,pz])
  (c,normal) = md_math.fitplane(points)

  # Create virtual atom in direction of normal
  v1 = Atom()
  v1.typenum = 0
  v1.type    = a[id].type + cnttypes
  v1.charge  = -0.100
  v1.x       = c[0] + 0.65*normal[0] # cg1-cge bond length = 0.6500
  v1.y       = c[1] + 0.65*normal[1] 
  v1.z       = c[2] + 0.65*normal[2]
  a[vaid] = v1 # Add atom to atom dictionary
  
  # Create new bond for virtual atom
  b1 = Bond()
  b1.type = 2
  b1.atomids = [id,vaid]
  b[vbid] = b1
  
  vaid += 1
  vbid += 1

  # Create virtual atom in opposite direction
  v2 = Atom()
  v2.typenum = 0
  v2.type    = a[id].type + cnttypes
  v2.charge  = -0.100
  v2.x       = c[0] - 0.65*normal[0] # cg1-cge bond length = 0.6500
  v2.y       = c[1] - 0.65*normal[1] 
  v2.z       = c[2] - 0.65*normal[2]
  a[vaid] = v2 # Add atom to atom dictionary

  # Create new bond for virtual atom
  b2 = Bond()
  b2.type = 2
  b2.atomids = [id,vaid]
  b[vbid] = b2

  vaid += 1
  vbid += 1


####################################
### Determine Angles & Dihedrals ###
####################################

# Find all angles & dihedrals
m.angles    = locality.findangles(m.bonds)
m.dihedrals = locality.finddihedrals(m.bonds,m.angles) 

an = m.angles

# Assign correct type for each angle
for id in an:
  atoms = an[id].atomids
  cntfs = []

  for i in range(0,3):
    if a[atoms[i]].type <= cnttypes:
      cntfs.append(1)      
    else:
      cntfs.append(0)

  if cntfs == [1,1,1]:
    an[id].type = 1
  elif cntfs == [0,1,1]:
    an[id].type = 2
  elif cntfs == [1,1,0]:
    an[id].type = 2
  elif cntfs == [0,1,0]:
    an[id].type = 3
  else:
    raise Exception("Angle type cannot be assigned")


# Count number of atoms,bonds,angles,dihedrals,impropers
m.natoms         = len(m.atoms)
m.natomtypes     = len(m.masses)
m.nbonds         = len(m.bonds)
m.nbondtypes     = len(m.bond_coeffs)
m.nangles        = len(m.angles)
m.nangletypes    = len(m.angle_coeffs)
m.ndihedrals     = len(m.dihedrals)
m.ndihedraltypes = len(m.dihedral_coeffs)
m.nimpropers     = 0
m.nimpropertypes = len(m.improper_coeffs)

write_pcffiff.moleculefile(file_out,m) 
