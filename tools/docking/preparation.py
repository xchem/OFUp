from ..misc import obconv
from htmd.ui import *
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers
from random import randint
import pybel


def prep_ligands(lig_sdf):
    '''
    convert ligands in an sdf file to pdbqt format
    :param lig_sdf: sdf file containing ligands
    :return:
    '''
    obconv(in_type='sdf', out_type='pdbqt', in_file=lig_sdf, out_file=lig_sdf.replace('.sdf', '.pdbqt'), options=['h'])

    return lig_sdf.replace('.sdf', '.pdbqt')


def prep_and_split(lig_sdf):
    suppl = Chem.SDMolSupplier(lig_sdf)
    no_mols = len([x for x in suppl])
    i = 1
    for mol in pybel.readfile("sdf", lig_sdf):
        mol.write("pdbqt", "%s_prepared.pdbqt" % lig_sdf.replace('.sdf', str(i)))
        i += 1

    names = [lig_sdf.split('.')[0] + str(n) + '_prepared.pdbqt' for n in range(1, no_mols + 1)]

    return names


def prep_protein(protein_pdb):
    '''
    use htmd and openbabel to prepare protein. Removes any non-protein atoms/mols
    :param protein_pdb: the path to the protein pdb file to use
    :return: prot: the htmd protein object
             preppdb: path to the prepared pdbqt file
    '''
    mol = Molecule(protein_pdb)
    mol.filter('protein')
    prot = proteinPrepare(mol)
    preppdb = protein_pdb.replace('.pdb', '_prepared.pdb')
    prot.write(preppdb)
    obconv(in_type='pdb', out_type='pdbqt', in_file=preppdb, out_file=preppdb.replace('.pdb', '.pdbqt'),
           options=['x', 'r'])

    return prot, preppdb.replace('.pdb', '.pdbqt')


def define_vina_box(ref_lig, output, bufx=10, bufy=10, bufz=10, exhaustiveness=20, seed=None, cpu=4):
    '''
    Use a reference ligand to produce a box definition for autodock vina
    :param ref_lig: a mol or sdf file of the reference ligand
    :param output: the name of the output file
    :param bufx: distance (A) of buffer in x direction (default=10)
    :param bufy: distance (A) of buffer in y direction (default=10)
    :param bufz: distance (A) of buffer in z direction (default=10)
    :param exhaustiveness: exhaustivness of the vina calculation (default=20)
    :param seed: random seed for vina (default=None)
    :param cpu: number of cpu's to use in vina calculation
    :return:
    '''
    mol = Chem.MolFromMolFile(ref_lig)
    conf = mol.GetConformer()
    params = rdShapeHelpers.ComputeConfBox(conf)

    c1 = np.array(params[0])
    c2 = np.array(params[1])

    center = np.mean((c1, c2), axis=0)

    dims = np.abs(c1 - c2)

    box_dims = [dims[0] + bufx, dims[1] + bufy, dims[2] + bufz]

    if not seed:
        seed = randint(0, 2147483647)

    with open(output, 'w') as f:
        f.write(
            """size_x = {}
    size_y = {}
    size_z = {}
    center_x = {}
    center_y = {}
    center_z = {}
    num_modes = 9999
    energy_range = 9999
    exhaustiveness = {}
    cpu={}
    seed={}
    """.format(box_dims[0], box_dims[1], box_dims[2],
               center[0], center[1], center[2],
               exhaustiveness,
               cpu,
               seed)
        )













