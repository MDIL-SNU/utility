"""
Code for converting VASP structure file to LAMMPS structure file

Usage:

python VASP_to_lammps.py [atom_type] [VASP structure file(POSCAR)] [LAMMPS structure file(somefilename.dat)]

* atom_type:
   ex) Si,O

ase, numpy is required to use this code.

Created by Kyuhyun Lee (190514)

Currently, MD trajectory is not supported.
"""

from ase import io
import numpy as np
import argparse

def rotAxis(cell, rot):
    if abs(cell[rot[0], rot[1]]) > 1e-4:
        vec1 = [cell[rot[0], rot[0]], cell[rot[0], rot[1]]]
        vec2 = [np.linalg.norm(vec1), 0.]
        cval = np.dot(vec1, vec2)/np.linalg.norm(vec1)/np.linalg.norm(vec2)
        sval = np.sqrt(1-cval**2)

        res_cell = np.copy(cell)
        for i in range(3):
            if cell[rot[0], rot[1]] > 0.:
                tmp_vec = [ res_cell[i,rot[0]]*cval + res_cell[i,rot[1]]*sval, \
                           -res_cell[i,rot[0]]*sval + res_cell[i,rot[1]]*cval]
            else:
                tmp_vec = [ res_cell[i,rot[0]]*cval - res_cell[i,rot[1]]*sval, \
                            res_cell[i,rot[0]]*sval + res_cell[i,rot[1]]*cval]

            res_cell[i,rot[0]] = tmp_vec[0]
            res_cell[i,rot[1]] = tmp_vec[1]

        return res_cell
    else:
        return cell

parser = argparse.ArgumentParser()
parser.add_argument('atom_types', \
                    help="Atom type seperated by ','. ex) 'Si,O' or 'Si,O,H' etc")
parser.add_argument('infile', \
                    help="VASP POSCAR/CONTCAR file.")
parser.add_argument('outfile', \
                    help="LAMMPS structure file.")
parser.add_argument('-f', '--format', default='vasp',
                    help="ASE io format string")
argv = parser.parse_args()

atoms = io.read(argv.infile, format=argv.format)
OUT = open(argv.outfile, 'w')

cell = atoms.get_cell()
scaled_posi = atoms.get_scaled_positions()

cell = rotAxis(cell, [0,1])
cell = rotAxis(cell, [0,2])
cell = rotAxis(cell, [1,2])

reposi = np.dot(scaled_posi, cell)

atom_types = argv.atom_types.split(',')
at_dict = dict()
masses = atoms.get_masses()
symbols = atoms.get_chemical_symbols()
m_dict = dict()

for i,item in enumerate(atom_types):
    at_dict[item] = i+1

    for j in range(len(symbols)):
        if symbols[j] == item:
            m_dict[item] = masses[j]
            break

print at_dict
print m_dict

OUT.write('LAMMPS structure file\n\n')
OUT.write(str(len(reposi)) + ' atoms\n' + \
          str(len(atom_types)) + ' atom types\n\n')
OUT.write('0.0 ' + str(cell[0,0]) + ' xlo xhi\n' + \
          '0.0 ' + str(cell[1,1]) + ' ylo yhi\n' + \
          '0.0 ' + str(cell[2,2]) + ' zlo zhi\n' + \
          str(cell[1,0]) + ' ' + str(cell[2,0]) + ' ' + str(cell[2,1]) + ' xy xz yz\n\n')
OUT.write('Masses\n\n')

for item in atom_types:
    OUT.write(str(at_dict[item]) + ' ' + str(m_dict[item]) + '\n')
OUT.write('\nAtoms\n\n')

for i,item in enumerate(reposi):
    OUT.write(str(i+1) + ' ' + str(at_dict[symbols[i]]) + ' ' + \
              ' '.join(item.astype(np.str)) + ' 0 0 0\n')

print cell
OUT.close()
