from rdkit import Chem


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