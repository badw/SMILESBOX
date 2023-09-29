from ase.build import molecule
import ase.data.pubchem as pubchem
from ase import Atoms
from ase.io import read,write

from tqdm import tqdm
import warnings
from shutil import copyfile
import os,glob,wget
from copy import deepcopy

class SMILESbox:
    
    def __init__(self):
        '''class for downloading and creating simulation-ready files for both xyz, .vasp, .cif'''
    
    def create_simulation_box(self,molecule=None,boxsize=10):
        mol = deepcopy(molecule)
        mol.set_cell([boxsize,boxsize,boxsize])
        mol.center()
        return(mol)
    
    def save_simulation_box(self,simulation_box=None,path=None,name=None):
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
    
    def save_xyz(self,molecule=None,path=None,name=None):
        if not name == None:
            if not path == None:
                end = os.path.join(str(path),str(name))
                molecule.write(filename=end)
            else:
                molecule.write(filename=name)
        else:
            if not path==None:
                end = os.path.join(str(path),str('molecule.xyz'))
                molecule.write(filename=end)
            else:
                molecule.write(filename='molecule.xyz') 
                
    def download_molecule_from_smiles(self,molecules=None):
        
        def _download_xyz(smiles):
            url = 'https://cactus.nci.nih.gov/chemical/structure/{}/file?format=xyz'.format(smiles)
            if not glob.glob('temp*.xyz') == []: # tidy up 
                [os.remove(x) for x in glob.glob('temp*.xyz')]
            wget.download(url,out='temp.xyz')
            molecule = read('temp.xyz')
            if not glob.glob('temp*.xyz') == []: # tidy up 
                [os.remove(x) for x in glob.glob('temp*.xyz')]
            return(molecule)
        
        downloaded_data = {}
        for mol,smiles in molecules.items():
            try:
                downloaded_data[mol] = _download_xyz(smiles)
            except:
                downloaded_data[mol] = None
        return(downloaded_data)