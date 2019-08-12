from ProLIF.prolif.ligand import Ligand
from ProLIF.prolif.protein import Protein
from ProLIF.prolif.fingerprint import Fingerprint
from rdkit import Chem
import rdkit.Chem.PandasTools as pdt


class PlifScore:
    def __init__(self, protein_pdb, reference_file):
        self.protein = Protein(inputFile=protein_pdb)
        self.protein.residuesFromMOL2File()

        suppl = Chem.SDMolSupplier(reference_file)
        mols = [Chem.AddHs(x) for x in suppl]
        ref = mols[0]
        self.reference_ligand = Ligand(ref)

        Fingerprint().generateIFP(ligand=self.reference_ligand, protein=self.protein)

    def score_conformers(self, file, write=False, method='tanimoto'):
        df = pdt.LoadSDF(file, embedProps=True)

        scores = []

        for i in range(0, len(df)):

            m = df.iloc[i]['ROMol']

            Chem.AddHs(m)
            ligand = Ligand(m)
            try:
                Fingerprint().generateIFP(ligand=ligand, protein=self.protein)
                ligand_score = round(float(ligand.getSimilarity(reference=self.reference_ligand, method=method)), 2)
                scores.append(ligand_score)
            except:
                scores.append(None)
            # df.iloc[i]['PLIF_SCORE'] = str(ligand_score)

        if write:
            pdt.WriteSDF(df, file)

        df['PLIF_SCORE'] = scores


        return df



