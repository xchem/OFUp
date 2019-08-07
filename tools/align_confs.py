from rdkit import Chem


def align_to_ref(to_align_file, hit_file):
    '''
    Align multiple conformers (or molecules) in a single sdf file to a reference mol in an sdf file
    :param to_align_file: the sdf file containing molecules to align
    :param hit_file: the sdf file containing the reference mol (uses position 0)
    :return:
    '''
    suppl = Chem.SDMolSupplier(hit_file)
    hit_mol = suppl[0]

    suppl = Chem.SDMolSupplier(to_align_file)
    mols = [x for x in suppl]

    for mol in mols:
        o3d = Chem.rdMolAlign.GetO3A(prbMol=mol, refMol=hit_mol)
        o3d.Align()

