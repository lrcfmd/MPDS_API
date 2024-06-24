import ase.io.vasp
import ase.io
import sys

filename=sys.argv[1]
atoms=ase.io.read(filename)
#ase.io.write('input.cif', atoms)
#ase.io.vasp.write_vasp('POSCAR', atoms, direct=True)
#atoms.write('POSCAR', format='vasp')
