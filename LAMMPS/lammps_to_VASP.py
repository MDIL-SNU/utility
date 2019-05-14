"""
Code for converting LAMMPS trajectory file to VASP XDATCAR and OUTCAR

Usage:

python lammps_to_VASP.py [atom_type] [LAMMPS trajectory(dump) file(somefilename.lammpstrj)]

* atom_type:
   ex) Si,O

ase, numpy is required to use this code.

Created by Kyuhyun Lee (190514)

This code converts LAMMPS dump file to VASP XDATCAR and OUTCAR file.
VASP OUTCAR file is not complete form.
It only contains atomic coordinates and atomic force.
(part starting from 'TOTAL-FORCE' in VASP original OUTCAR)

If your structure contains 4 or more atom types, 
add atom type in 'dummy_sym' list in the order of atomic number. 
"""

from ase import io
import numpy as np
import argparse

global dummy_sym
dummy_sym = ['H', 'He', 'Li']

def write_header(XDAT, atoms, atom_types):
    cells = atoms.get_cell()
    atom_nums = np.zeros(len(atom_types)).astype(np.int)
    XDAT.write(' {:12.10f} {:12.10f} {:12.10f}\n {:12.10f} {:12.10f} {:12.10f}\n {:12.10f} {:12.10f} {:12.10f}\n'.\
                format(cells[0,0], cells[0,1], cells[0,2], \
                       cells[1,0], cells[1,1], cells[1,2], \
                       cells[2,0], cells[2,1], cells[2,2]))

    symbols = atoms.get_chemical_symbols()
    #dummy_sym = ['H', 'He', 'Li']
    for item in symbols:
        for j,jtem in enumerate(atom_types):
            if item==dummy_sym[j]:
                atom_nums[j] += 1
                break

    atom_totnum = sum(atom_nums)

    XDAT.write('  ' + '   '.join(atom_types) + '\n')
    XDAT.write('  ' + '   '.join(atom_nums.astype(np.str)) + '\n')

    return atom_nums, atom_totnum

def write_posis(XDAT, atoms, atom_idx):
    tails = '\n'

    scaled_posi = atoms.get_scaled_positions()
    scaled_posi = scaled_posi[atom_idx, :]

    for item in scaled_posi:
        XDAT.write('  ' + ' '.join(item.astype(np.str)) + tails)

parser = argparse.ArgumentParser()
parser.add_argument('atom_types', \
                    help="Atom type seperated by ','. ex) 'Si,O' or 'Si,O,H' etc")
parser.add_argument('infile', \
                    help="lammps trajectory file.")
argv = parser.parse_args()

XDAT = open('XDATCAR', 'w')
OUT = open('OUTCAR', 'w')

atom_types = argv.atom_types.split(',')

XDAT.write('lammps to VASP\n1.0\n')
    
atoms = io.read(argv.infile, index=':', format='lammps-dump')
atom_nums, atom_totnum = write_header(XDAT, atoms[0], atom_types)
header_num = 9

for item in atom_types:
    OUT.write('POTCAR: PBE {}\n'.format(item))
for item in atom_types:
    OUT.write('POTCAR: PBE {}\n'.format(item))

tmp_str = '   ions per type =           '
for item in atom_nums:
    tmp_str += '  ' + str(item)

OUT.write(tmp_str + '\n')

with open(argv.infile, 'r') as fil:
    for i,item in enumerate(atoms):
        E_sum = 0.
        F_list = np.zeros([atom_totnum, 6])
        syms = np.array(item.get_chemical_symbols())
        atom_idx = np.where(syms == dummy_sym[0])[0]
        for jtem in dummy_sym[1:]:
            atom_idx = np.concatenate((atom_idx, np.where(syms == jtem)[0]), axis=0)

        for j in range(header_num):
            tmp_line = fil.readline()

        for j in range(atom_totnum):
            tmp_line = fil.readline().replace('\n','').split()
            E_sum += float(tmp_line[-1])
            #atom_idx[int(tmp_line[0])-1] = atom_nums[tmp_idx_type] + idx_pertype[tmp_idx_type]
            F_list[atom_idx[int(tmp_line[0])-1], :3] = np.array(tmp_line[2:5]).astype(np.float64)
            F_list[atom_idx[int(tmp_line[0])-1], 3:] = np.array(tmp_line[8:11]).astype(np.float64)

        if i != 0:
            XDAT.write('lammps to VASP\n1.0\n')
            _, _ = write_header(XDAT, item, atom_types)
        XDAT.write('Direct configuration= {:5}\n'.format(i+1))
        write_posis(XDAT, item, atom_idx)

        cell_info = item.get_cell()
        OUT.write('Iteration {:6}\n'.format(i+1))
        OUT.write('      direct lattice vectors\n')
        OUT.write('  {:12.8f} {:12.8f} {:12.8f}\n'.format(cell_info[0][0], cell_info[0][1], cell_info[0][2]))
        OUT.write('  {:12.8f} {:12.8f} {:12.8f}\n'.format(cell_info[1][0], cell_info[1][1], cell_info[1][2]))
        OUT.write('  {:12.8f} {:12.8f} {:12.8f}\n'.format(cell_info[2][0], cell_info[2][1], cell_info[2][2]))
        OUT.write(' POSITION                                       TOTAL-FORCE (eV/Angst)\n')
        OUT.write(' -----------------------------------------------------------------------------------\n')
        for jtem in F_list:
            OUT.write(' {:12.5f} {:12.5f} {:12.5f}   {:13.6f} {:13.6f} {:13.6f}\n'.\
                format(jtem[0], jtem[1], jtem[2], jtem[3], jtem[4], jtem[5]))
        OUT.write(' -----------------------------------------------------------------------------------\n')
        OUT.write('  free  energy   TOTEN  = {:18.8f} eV\n\n\n'.format(E_sum))

XDAT.close()
OUT.close()
