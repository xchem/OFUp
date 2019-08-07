from rdkit import Chem
import os
from rdkit.Chem import AllChem, rdShapeHelpers


def generateconformations(m, n):
    '''
    Generate n conformations of rdkit mol (m) using ETKDG method
    :param m: the rdkit molecule object to generate conformers for
    :param n: the number of conformations to generate and embed into m
    :return: the mol with the conformations embedded; list of ids of the conformers
    '''
    # add hydrogen to molecule
    m = Chem.AddHs(m)
    # generate multiple conformations and return list of conf ids
    ids = AllChem.EmbedMultipleConfs(m, numConfs=n, params=AllChem.ETKDG())
    # align all of the conformers to each other
    _ = Chem.rdMolAlign.AlignMolConformers(m)
    # return the mol object with embedded conformers, and list of their ids
    return m, list(ids)


def writeout_confs(outname, mol, ids):
    '''
    Write conformers from an rdkit mol object out into an sdf file.
    :param outname: the name of the file to write out to
    :param mol: the rdkit mol object
    :param ids: a list of ids for conformers to write out
    :return: returns the name of the file written
    '''
    writer = Chem.SDWriter(outname)
    for i in ids:
        writer.write(mol, confId=i)

    return outname


def confs_and_write(candidates_file, no_confs):
    '''
    Take an SDF file of molecules and generate n conformers for each molecule.
    The resulting structures will be written out to a file named <original-name>_<index-in-sdf>_confs.sdf

    NB: an output file will not be overwritten if it already exists.

    :param candidates_file: input sdf file of molecules
    :param no_confs: the number of conformations to generate for each molecule
    :return: a list of the files generated
    '''
    # read in all candidate mols
    suppl = Chem.SDMolSupplier(candidates_file)
    mols = [x for x in suppl]
    out_list = []
    # for every candidate
    for i in range(0, len(mols)):
        # name out file by index of mol from initial file
        out_file = candidates_file.replace('.sdf', str('_' + str(i) + '_confs.sdf'))
        # if this candidate hasn't been done
        if not os.path.isfile(out_file):
            print('processing ' + candidates_file + ': candidate ' + str(i) + '(total=' + str(len(mols)) + ')')
            # generate n=no_confs conformers and write to out_file
            mol, ids = generateconformations(mols[i], no_confs)
            o = writeout_confs(outname=out_file, mol=mol, ids=ids)
            out_list.append(o)
        # if done, move onto next molecule
        else:
            continue

    return out_list

