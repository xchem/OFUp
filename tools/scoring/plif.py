from .ProLIF.prolif.ligand import Ligand
from .ProLIF.prolif.protein import Protein
from .ProLIF.prolif.fingerprint import Fingerprint
from rdkit import Chem
import rdkit.Chem.PandasTools as pdt


class PlifScore:
    def __init__(self, protein_pdb, reference_file):
        self.protein = Protein(inputFile=protein_pdb)
        self.protein.residuesFromMOL2File()

        if reference_file:
            suppl = Chem.SDMolSupplier(reference_file)
            mols = [Chem.AddHs(x) for x in suppl]
            ref = mols[0]

            self.reference_ligand = Ligand(ref)

            Fingerprint().generateIFP(ligand=self.reference_ligand, protein=self.protein)

    def score_conformers(self, file, write=False, method='tanimoto', score_col='PLIF_SCORE'):
        df = pdt.LoadSDF(file, embedProps=True)

        scores = []

        for i in range(0, len(df)):

            m = df.iloc[i]['ROMol']

            Chem.AddHs(m)
            ligand = Ligand(m)
            try:
                Fingerprint().generateIFP(ligand=ligand, protein=self.protein)
                ligand_score = round(float(ligand.getSimilarity(reference=self.reference_ligand, method=method)), 2)
                scores.append(str(ligand_score))
            except:
                scores.append(None)

        df[score_col] = scores

        if write:
            pdt.WriteSDF(df, file, properties=list(df.columns))

        return df

    def score_mol(self, mol):
        self.lig = Ligand(mol)
        self.lig.ifp = Fingerprint().generateIFP(ligand=self.lig, protein=self.protein)

        return self.lig.ifp


def cluster_vects(fps, cutoff=0.2):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


def cluster_plifs(mol_list, protein, cutoff=0.5):

    prot_plif = PlifScore(protein_pdb=protein, reference_file=None)
    vects = [prot_plif.score_mol(lig) for lig in mol_list if lig]
    clusters = cluster_vects(vects, cutoff=cutoff)

    return clusters






