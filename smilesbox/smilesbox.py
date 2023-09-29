from ase.build import molecule
import ase.data.pubchem as pubchem
from ase import Atoms
from ase.io import read,write

from tqdm import tqdm
import warnings
from shutil import copyfile
import os,glob,wget
from copy import deepcopy
from openbabel import pybel

class SMILESbox:
    
    def __init__(self):
        '''class for downloading and creating simulation-ready files for both xyz, .vasp, .cif'''
    
    def create_simulation_box(self,molecule=None,boxsize=10):
        mol = deepcopy(molecule)
        mol.set_cell([boxsize,boxsize,boxsize])
        mol.center()
        return(mol)
    
    def rotate(self,which='x',amount=90,reset=True):

    
    def save_vasp(self,simulation_box=None,path=None,name=None):
        if not name == None:
            if not path == None:
                end = os.path.join(str(path),str(name))
                simulation_box.write(filename=end,vasp5=True,sort=True)
            else:
                simulation_box.write(filename=name,vasp5=True,sort=True)
        else:
            if not path==None:
                end = os.path.join(str(path),str('POSCAR'))
                simulation_box.write(filename=end,vasp5=True,sort=True)
            else:
                simulation_box.write(filename='POSCAR',vasp5=True,sort=True)             
    
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
        self.atoms == atoms
        #return(atoms)
