import numpy
import pymatgen
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element

Rneib=3.6
V=63.81686774
Conv=16.021766

#Read and analyse structure provided
struct = pymatgen.Structure.from_str(open("bulk_tetragonal_CONTCAR.vasp").read(), fmt="poscar")

#Calculate dipole moment and polarisation for each Ti in the structure
#
for site in struct.sites: #First look over all sites in structure
  if(site.specie==Element("Ti")): #If site is a Ti
    p=4.0*site.coords #Calculate dipole moment contribution from Ti
    neigbours=struct.get_neighbors(site,Rneib) #find neigbouring sites within Rneib of Ti
    for neib in neigbours: #For each of the neighbouring sites
      if(neib[0].specie==Element("Ba")): #If it is a Ba add the dipole moment contribution from Ba
        p=p+0.25*neib[0].coords 
      if(neib[0].specie==Element("O")):  #If it is a O add the dipole moment contribution from O
        p=p-1.0*neib[0].coords
    P=Conv*p/V #Convert to P in units of C/m2
    print(site.coords[0],P) #Print x-coordinate of central Ti and the associated polarisation vector
