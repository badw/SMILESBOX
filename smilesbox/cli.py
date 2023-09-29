import argparse
import os
import sys
from smilesbox.smilesbox import SMILESbox
import logging
import numpy as np 

__author__ = "Benjamin A. D. Williamson"
__version__ = "1.0"
__maintainer__ = "Benjamin A. D. Williamson"
__email__ = "benjamin.williamson@ntnu.no"
__date__ = "September 2023"

def smilesbox(
        filename="POSCAR",
        directory=".",
        save=False,
        smiles=None,
        box=True,
        dimensions=[10,10,10],
        rotate=False,
        rotate_vector='x',
        rotate_angle='90'
        ):
    sb = SMILESbox()
    if smiles==None:
         logging.error(
              'no SMILES given :('
              )
    else:
        sb.smiles_to_atoms(smiles)
        
    if rotate == True:
        sb.rotate(v=rotate_vector,a=rotate_angle)
    if box == True:
        sb.add_box(dimensions=dimensions)
    if save == True:
        sb.save(filename=filename,directory=directory)

    logging.info(
             
             '''SMILES:

 {}

FORMULA:

 {}

CELL:

{}

{}'''.format(smiles,
                                   sb.atoms.get_chemical_formula(),
                                   str(np.array(sb.atoms.cell)),
                                   sb.molecule.write(format='ascii',opt={'w':50,'h':20,'m':0,'a':2})
                                   )
             )

def _get_parser():
    parser = argparse.ArgumentParser(
        description="""smilesbox generates a simulation cell (either .xyz or .vasp) from a SMILES string""",
        epilog=f"""Author: {__author__}
    Version:{__version__}
    Last updated: {__date__}""",
    )
        
    parser.add_argument(
        "-s",
        "--smiles",
        default=None,
        help="SMILES string i.e. HH(C)"
    ) 

    parser.add_argument(
        "--save",
        action="store_true",
        help="save to file? default=False",
    )

    parser.add_argument(
        "-f",
        "--filename",
        default="POSCAR",
        help="chosen filename"
    )

    parser.add_argument(
        "-d",
        "--directory",
        default=".",
        help="directory in which to save (default='.')"
    )

    parser.add_argument(
        "-r",
        "--rotate",
        action="store_true",
        help="rotate? default=False"
    )

    parser.add_argument(
        "-v",
        "--vector",
        default="x",
        help="rotation vector: x,y, or z"
    )

    parser.add_argument(
        "-a",
        "--angle",
        default=90,
        help="rotation angle, default=90 degrees"
    )

    parser.add_argument(
        "-b",
        "--box",
        action="store_true",
        help="add simulation box? default=True"
    )

    parser.add_argument(
        "--dimensions",
        default=[10,10,10],
        help="simulation box dimensions. default=[10,10,10]"
    )
    return(parser)

def main():
    args = _get_parser().parse_args()

    logging.basicConfig(
        filename="smilesbox.log",
        level=logging.DEBUG,
        filemode="w",
        format="%(message)s"
    )
    console=logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger("").addHandler(console)

    smilesbox(
        filename=args.filename,
        directory=args.directory,
        save=args.save,
        smiles=args.smiles,
        box=args.box,
        dimensions=args.dimensions,
        rotate=args.rotate,
        rotate_vector=args.vector,
        rotate_angle=args.angle
)
    
if __name__ == "__main__":
    main()
    
