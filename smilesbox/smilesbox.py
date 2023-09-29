from ase.build import molecule
import ase.data.pubchem as pubchem
from ase import Atoms
import ase
from ase.io import read,write

from tqdm import tqdm
import warnings
from shutil import copyfile
import os,glob,wget
from copy import deepcopy
from openbabel import pybel
import numpy as np 

class SMILESbox:
    
    def __init__(self):
        '''class for downloading and creating simulation-ready files for both xyz, .vasp, .cif'''
    
    def add_box(self,dimensions=[10,10,10]):
        '''default = [10,10,10]
        todo:
        * add whether molecule has enough vacuum'''
        self.atoms.set_cell(dimensions)
        self.atoms.set_pbc(True)
        self.atoms.center()
    
    def rotate(self,vector='x',angle=90):
        self.atoms.rotate(v=vector,angle=angle)
    
    def save(self,filename='POSCAR',directory='.'):
        ''' checks if there is a box installed if it prints to vasp'''
        if not np.sum(self.atoms.cell) == 0:
            self.atoms.write(os.path.join(directory,'{}{}'.format(filename,'.vasp')),vasp5=True,sort=True)        
        else:
            self.atoms.write(os.path.join(directory,'{}{}'.format(filename,'.xyz')),xyz=True)  
    
    def smiles_to_atoms(self,smiles=None):
        molecule = pybel.readstring("smi",smiles)
        molecule.make3D()
        symbols = []
        positions = []
        for i in molecule.atoms:
            coords = list(i.coords)
            atom = [k for k,v in ase.data.atomic_numbers.items() if v == i.atomicnum]
            symbols.append(atom[0])
            positions.append(coords)
        atoms = Atoms(symbols=symbols,positions=positions)
        self.atoms = atoms
        self.molecule=molecule
        #return(atoms)
