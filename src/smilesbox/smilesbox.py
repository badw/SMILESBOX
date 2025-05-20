from ast import Attribute
from ase import Atoms
from ase.data import atomic_numbers
from ase.io import read,write
import os
from openbabel import pybel
import numpy as np 
import math 

class SMILESbox:
    
    def __init__(self,smiles:str):
        
        self.smiles = smiles 
        self.atoms = self._smiles_to_atoms()


    def _smiles_to_atoms(
            self,
            )->Atoms:
        
        molecule = pybel.readstring("smi",self.smiles)
        molecule.make3D()
        symbols = []
        positions = []
        for i in molecule.atoms:
            coords = list(i.coords)
            atom = [k for k,v in atomic_numbers.items() if v == i.atomicnum]
            symbols.append(atom[0])
            positions.append(coords)

        atoms = Atoms(symbols=symbols,positions=positions)
        self.molecule = molecule
        return(atoms)
    

    def generate_random_points(self,N, L, min_dist, max_attempts=100000):
        """
        Generate N random 3D points inside a cube of size L, with a minimum spacing.    

        Parameters:
        - N: Number of points.
        - L: Length of the cube along each axis (cube from 0 to L in x, y, z).
        - min_dist: Minimum allowed distance between any two points.
        - max_attempts: Max number of iterations before giving up.    

        Returns:
        - points: (N, 3) array of valid 3D coordinates.
        """
        points = []
        attempts = 0    

        while len(points) < N and attempts < max_attempts:
            candidate = np.random.uniform(0, L, size=3)
            if all(np.linalg.norm(candidate - np.array(p)) >= min_dist for p in points):
                points.append(candidate)
            attempts += 1    

        if len(points) < N:
            raise RuntimeError(f"Could not place all points within {max_attempts}     attempts. "
                               f"Try reducing N or min_dist.")
        return np.array(points)
    
    def multiple_in_a_box(self,num_points,box_size,min_distance=0.5,random_rotate=True):
        final_cell = Atoms(cell=[box_size,box_size,box_size])
        try:
            points = self.generate_random_points(N=num_points,L=box_size,min_dist=min_distance)
        except RuntimeError as e:
            raise RuntimeError(f'{e} or try a larger box...{box_size} may be too small')
            
        for point in points:
            _atoms = self._smiles_to_atoms()
            _atoms = self.add_box(_atoms,[box_size,box_size,box_size])
            if random_rotate:
                _atoms = self.rotate(_atoms,np.random.choice(['x','y','z']),np.random.randint(360))
            _atoms = self.translate(_atoms,point)
            final_cell+= _atoms
        return(final_cell)
    
    def add_box(
            self,
            atoms:Atoms=None,
            dimensions:np.ndarray=None,
            )->None:
        """
        dimensions = np.ndarray([[10,0,0],
        [0,10,0],
        [0,0,10]])

        i.e. 
        import numpy as np 
        dimensions = np.eye(3) * 10 
        """
        if not atoms:
            atoms = self.atoms    
        atoms.set_cell(dimensions)
        atoms.set_pbc(True)
        atoms.center()
        self.atoms = atoms
        return(atoms)


    
    def translate(self,atoms,new_coords):
        com = atoms.get_center_of_mass()
        vector = new_coords - com 
        atoms.translate(vector)
        return(atoms)

    
    def rotate(self,atoms,vector='x',angle=90):
        atoms.rotate(a=angle,v=vector)
        return(atoms)


