import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit import RDConfig

# SuCOS code adapted from:
# https://bitbucket.org/Susanhleung/sucos/src/master/calc_SuCOS.py


class SuCOS:
    def __init__(self, reference=None, target_s=None):
        '''
        Setup of the SuCOS scoring class. If you do not set reference or target_s, you will have to feed them
        as rdkit mol objects to each class method/function.

        Sets: self.reference: an rdkit mol object of the reference mol
              self.target_s: a list of rdkit mol objects for the target mols

        :param reference: an sdf or mol file containing the reference (hit) molecule
        :param target_s: an sdf or mol file containing one or multiple molecules to compare against reference
        '''
        self.fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

        self.fmParams = {}
        for k in self.fdef.GetFeatureFamilies():
            self.fparams = FeatMaps.FeatMapParams()
            self.fmParams[k] = self.fparams

        self.keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
                'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')

        if reference:
            refs = Chem.SDMolSupplier(reference, sanitize=True)
            self.reference = [x for x in refs if x][0]

        if target_s:
            tgts = Chem.SDMolSupplier(target_s, sanitize=True)
            self.target_s = [x for x in tgts if x]


    def get_fm_score(self, small_m=None, large_m=None,
                            score_mode=FeatMaps.FeatMapScoreMode.All):
        '''
        Get the feature-map score for a small molecule vs. a large molecule
        :param small_m: rdkit mol object for the small molecule
        :param large_m: rdkit mol object for the large molecule
        :param score_mode: default = featuremaps score, defined in init
        :return: the value of the feature map score
        '''

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
        '''
        Get the SuCOS score for one mol compared to another mol (mol=rdkit mol object)
        :param mol1: rdkit mol object (small mol)
        :param mol2: rdkit mol object (large mol)
        :param score_mode: default = featuremaps score, defined in init
        :return:
        '''
        ref = Chem.AddHs(mol1)
        prb = Chem.AddHs(mol2)

        fm_score = self.get_fm_score(ref, prb, score_mode)
        fm_score = np.clip(fm_score, 0, 1)

        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref, prb,
                                                         allowReordering=False)
        protrude_dist = np.clip(protrude_dist, 0, 1)

        SuCOS_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

        return SuCOS_score


    def sucos_mol_to_many(self, ref_file, prb_file, score_mode=FeatMaps.FeatMapScoreMode.All, p=False):
        '''
        compare multiple mol objects to one reference mol
        :param ref_file: sdf or mol file of reference (hit) mol
        :param prb_file: sdf file of mols to compare to reference
        :param score_mode: default = featuremaps score, defined in init
        :param p: True/False - print sucos and tanimoto score for every comparison
        :return:
            results_sucos: a dict containing idx of each prb mol as key, and sucos score as value
            results_tani: a dict containing idx of each prb mol as key, and tanimoto similarity distance score as value
            smi_mol: a list of smiles strings for the prb mols
            prb_mols: a list of rdkit mol objects for the prb mols
            reflig: an rdkit mol object for the reference
        '''

        reflig = Chem.MolFromMolFile(ref_file, sanitize=True)
        prb_mols = Chem.SDMolSupplier(prb_file, sanitize=True)
        prb_mols = [x for x in prb_mols if x]
        idx = 0

        results_sucos = {}
        results_tani = {}

        smi_mol = Chem.MolToSmiles(prb_mols[0])

        for prb_mol in prb_mols:

            SuCOS_score = self.sucos_mol_to_mol(mol1=reflig, mol2=prb_mol, score_mode=score_mode)
            tanimoto_score = Chem.rdShapeHelpers.ShapeTanimotoDist(reflig, prb_mol)

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
