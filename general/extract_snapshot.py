"""
Code for extracting snapshot from VASP or LAMMPS trajectory file.

Usage:

python extract_snapshot.py [trajectory_filename] [index] [format]

* index:
   start:end:interval
   ex) 1:10:2 -> 1,3,5,7,9
       :10 -> 1,2,3,4,5,6,7,8,9
       : -> all
* format:
   ase file format.
   vasp: VASP POSCAR
   vasp-out: VASP OUTCAR
   lammps-dump: LAMMPS trajectory file

ase is required to use this code.

Created by Kyuhyun Lee (190514)
"""

import sys
from ase import io

atoms = io.read(sys.argv[1], index=sys.argv[2], format=sys.argv[3])

for i,item in enumerate(atoms):
    io.write('POSCAR_T{}'.format(i+1), item, format='vasp')
