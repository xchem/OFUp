from preparation import prep_and_split, prep_protein
import shutil
import os
import subprocess
from joblib import Parallel, delayed
import multiprocessing

vina_exe = shutil.which('vina')

def vina_run(prot, lig, conf=None, cx=None, cy=None, cz=None, sx=None, sy=None, sz=None, exh=None, cpus=4, exe=vina_exe):
    if not vina_exe:
        raise Exception('vina executable was not found on your path, please specify it!')
    out = lig.replace('.sdf', '_docked.sdf')
    log = lig.replace('.sdf', '_docked.log')
    if not conf:
        command = f"{exe} --center_x {cx} --center_y {cy} --center_z {cz} --size_x {sx} --size_y {sy} --size_z {sz}" \
            f"--exhaustiveness {exh} --num_modes 9999 --energy_range 9999 --receptor {prot} --ligand {lig}" \
            f"--out {out} --log {log} --cpu {cpus}"
    else:
        command = f"{exe} --conf {conf} --num_modes 9999 --energy_range 9999 --receptor {prot} --ligand {lig}" \
            f"--out {out} --log {log}"

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()

    return out, err



def vina_dock_multi(protein_pdb, ligands_sdf,
                    conf=None, docking_dir=os.getcwd(), multiproc=False, exe=vina_exe):

    shutil.copy(ligands_sdf, docking_dir)
    shutil.copy(protein_pdb, docking_dir)

    new_ligs_path = os.path.join(docking_dir, ligands_sdf)
    new_prot_path = os.path.join(docking_dir, protein_pdb)

    lig_files = prep_and_split(new_ligs_path)
    _, protein_pdbqt = prep_protein(new_prot_path)

    if not multiproc:
        cwd = os.getcwd()
        os.chdir(docking_dir)
        for lig in lig_files:
            out, err = vina_run(prot=new_prot_path, lig=lig, conf=conf, exe=exe)
        os.chdir(cwd)

    else:
        cwd = os.getcwd()
        os.chdir(docking_dir)
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(vina_run)(prot=new_prot_path, lig=lig,
                                                       conf=conf, exe=exe) for lig in lig_files)
        os.chdir(cwd)