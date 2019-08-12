from rdkit import Chem
import openbabel
import pybel


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


def add_property(mol, prop_name, prop_value):
    newData = openbabel.OBPairData()
    newData.SetAttribute(prop_name)
    newData.SetValue(prop_value)
    mol.OBMol.CloneData(newData)


def vina_pdbqt_to_sdf(ifile, output):
    docked = pybel.readfile('pdbqt', ifile)
    sdf = pybel.Outputfile("sdf", output, overwrite=True)
    for mol in docked:
        if mol.OBMol.HasData('REMARK'):
            remark = mol.OBMol.GetData('REMARK').GetValue()
            lines = remark.splitlines()
            tokens = lines[0].split()

            # add the score property
            add_property(mol, "SCORE", tokens[2])
            # add the first RMSD property
            add_property(mol, "RMSD_LB", tokens[3])
            # add the second RMSD property
            add_property(mol, "RMSD_UB", tokens[4])

        sdf.write(mol)

    sdf.close()