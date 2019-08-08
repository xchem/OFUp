from rdkit import Chem
import openbabel


def write_mols(mols, outname):
    '''
    Write a list of rdkit mols out to an sdf with outname
    :param mols: list of rdkit mols
    :param outname: name of sdf file to write
    :return: name of file written
    '''

    writer = Chem.SDWriter(outname)
    for mol in mols:
        writer.write(mol)

    return outname


def obconv(in_type, out_type, in_file, out_file, options=[]):
    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats(in_type, out_type)
    for o in options:
        obc.AddOption(o)

    mol = openbabel.OBMol()
    obc.ReadFile(mol, in_file)
    mol.AddHydrogens()
    obc.WriteFile(mol, out_file)