from .preparation import prep_and_split, prep_protein
import shutil
import os
import subprocess
from joblib import Parallel, delayed
import multiprocessing
import glob
from tqdm import tqdm

vina_exe = shutil.which('vina')

def vina_run(prot, lig, conf=None, cx=None, cy=None, cz=None, sx=None, sy=None, sz=None, exh=None, cpus=4, exe=vina_exe):
    if not vina_exe:
        raise Exception('vina executable was not found on your path, please specify it!')
    out = lig.replace('.pdbqt', '_docked.pdbqt')
    log = lig.replace('.pdbqt', '_docked.log')
    if not conf:
        command = f"{exe} --center_x {cx} --center_y {cy} --center_z {cz} --size_x {sx} --size_y {sy} --size_z {sz}" \
            f"--exhaustiveness {exh} --num_modes 9999 --energy_range 9999 --receptor {prot} --ligand {lig}" \
            f"--out {out} --log {log} --cpu {cpus}"
    else:
        command = f"{exe} --config {conf} --receptor {prot} --ligand {lig}" \
            f" --out {out} --log {log}"

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()

    return out, err, command



def vina_dock_multi(protein_pdb, ligands_sdf,
                    conf=None, docking_dir=os.getcwd(), multiproc=False, exe=vina_exe):
    '''
    Dock multiple ligands from an sdf file against a protein, with a given configuration file with vina.
    Also does preparation from pdb and sdf

    :param protein_pdb: the input pdb file to prepare and dock against
    :param ligands_sdf: an sdf file of multiple ligands to split, prepare and dock
    :param conf: the path of a configuration file (required for now)
    :param docking_dir: the output directory
    :param multiproc: use multi-processing? (boolean)
    :param exe: the path to the vina executable
    :return:
    '''
    shutil.copy(ligands_sdf, docking_dir)
    shutil.copy(protein_pdb, docking_dir)
    if conf:
        shutil.copy(conf, docking_dir)
        conf = os.path.join(docking_dir, conf.split('/')[-1])
        os.chmod(conf, 775)

    new_ligs_path = os.path.join(docking_dir, ligands_sdf.split('/')[-1])
    os.chmod(new_ligs_path, 775)
    new_prot_path = os.path.join(docking_dir, protein_pdb.split('/')[-1])
    os.chmod(new_prot_path, 775)

    lig_files = prep_and_split(new_ligs_path)
    _, protein_pdbqt = prep_protein(new_prot_path)

    if not multiproc:
        cwd = os.getcwd()
        os.chdir(docking_dir)
        for lig in lig_files:
            out, err, command = vina_run(prot=protein_pdbqt, lig=lig, conf=conf, exe=exe)
        os.chdir(cwd)

    else:
        cwd = os.getcwd()
        os.chdir(docking_dir)
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(vina_run)(prot=protein_pdbqt, lig=lig,
                                                     conf=conf, exe=exe) for lig in tqdm(lig_files))
        os.chdir(cwd)
        
    if not os.path.isdir(os.path.join(docking_dir, 'docked/')):
        os.mkdir(os.path.join(docking_dir, 'docked/'))
    if not os.path.isdir(os.path.join(docking_dir, 'logs/')):
        os.mkdir(os.path.join(docking_dir, 'logs/'))
    if not os.path.isdir(os.path.join(docking_dir, 'inputs/')):
        os.mkdir(os.path.join(docking_dir, 'inputs/'))

    docked_files = glob.glob(os.path.join(docking_dir, '*_docked.pdbqt'))
    for f in docked_files:
        shutil.move(f, os.path.join(docking_dir, 'docked/'))

    logs = glob.glob(os.path.join(docking_dir, '*_docked.log'))
    for f in logs:
        shutil.move(f, os.path.join(docking_dir, 'logs/'))

    in_ligs = glob.glob(os.path.join(docking_dir, '*_prepared.pdbqt'))
    for f in in_ligs:
        shutil.move(f, os.path.join(docking_dir, 'in_ligs/'))
