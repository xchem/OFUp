import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit import RDConfig

# SuCOS code adapted from:
# https://bitbucket.org/Susanhleung/sucos/src/master/calc_SuCOS.py)


class SuCOS:
    def __init__(self, reference=None, target_s=None):
        self.fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

        self.fmParams = {}
        for k in self.fdef.GetFeatureFamilies():
            self.fparams = FeatMaps.FeatMapParams()
            self.fmParams[k] = self.fparams

        self.keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
                'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')

        if reference:
            self.reference = reference

        if target_s:
            self.target_s = target_s

    def get_fm_score(self, small_m=None, large_m=None,
                            score_mode=FeatMaps.FeatMapScoreMode.All):

        if not small_m and not self.reference:
            raise Exception('Reference molecule not set! Please specify or provide to init')

        if not large_m and not self.target_s:
            raise Exception('Target molecule(s) not set! Please specify or provide to init')

        if not small_m:
            small_m = self.reference

        if not large_m:
            large_m = self.target_s

        feat_lists = []

        for m in [small_m, large_m]:
            rawFeats = self.fdef.GetFeaturesForMol(m)
            # filter that list down to only include the ones we're interested in
            feat_lists.append([f for f in rawFeats if f.GetFamily() in self.keep])

        fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=self.fmParams) for x in feat_lists]
        fms[0].scoreMode = score_mode
        fm_score = fms[0].ScoreFeats(feat_lists[1]) / min(fms[0].GetNumFeatures(), len(feat_lists[1]))

        return fm_score

    def sucos_mol_to_mol(self, mol1, mol2, score_mode=FeatMaps.FeatMapScoreMode.All):
        ref = Chem.AddHs(mol1)
        prb = Chem.AddHs(mol2)

        fm_score = self.get_FeatureMapScore(ref, prb, score_mode)
        fm_score = np.clip(fm_score, 0, 1)

        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref, prb,
                                                         allowReordering=False)
        protrude_dist = np.clip(protrude_dist, 0, 1)

        SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

        return SuCOS_score

    def sucos_mol_to_many(self, ref_file, prb_file, score_mode=FeatMaps.FeatMapScoreMode.All, p=False):
        reflig = Chem.MolFromMolFile(ref_file, sanitize=True)
        ref = Chem.AddHs(reflig)
        prb_mols = Chem.SDMolSupplier(prb_file, sanitize=True)
        prb_mols = [x for x in prb_mols if x]
        idx = 0

        results_sucos = {}
        results_tani = {}

        smi_mol = Chem.MolToSmiles(prb_mols[0])

        for prb_mol in prb_mols:

            prb = Chem.AddHs(prb_mol)

            fm_score = get_FeatureMapScore(ref, prb, score_mode)
            fm_score = np.clip(fm_score, 0, 1)

            protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref, prb,
                                                             allowReordering=False)
            protrude_dist = np.clip(protrude_dist, 0, 1)

            SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)
            tanimoto_score = Chem.rdShapeHelpers.ShapeTanimotoDist(ref, prb)

            results_sucos[str(idx)] = SuCOS_score
            results_tani[str(idx)] = tanimoto_score

            if p:
                print ("********************************")
                print ("index: " + str(idx))
                print ("SuCOS score:\t%f" % SuCOS_score)
                print ("Tani score:\t%f" % tanimoto_score)
                print ("********************************")

            idx += 1

        return results_sucos, results_tani, smi_mol, prb_mols, reflig
