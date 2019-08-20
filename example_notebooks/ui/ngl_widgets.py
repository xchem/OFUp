from ipywidgets import BoundedIntText, Button, GridspecLayout, ToggleButton, Text
import nglview as nv
from rdkit import Chem
from nglview.show import StringIO


class interactive_view:

    def __init__(self, pdb_file, mol_df):
        self.pdb = Chem.MolFromPDBFile(pdb_file)
        self.view = nv.show_rdkit(self.pdb)
        self.view.stage.set_parameters(mouse_preset='coot')
        self.df = mol_df
        self.current_component = None

    def merge_mols(self, lig_mol):
        mol = Chem.CombineMols(self.pdb, lig_mol)
        return mol

    def add_mol(self, chg=None):
        if self.current_component:
            self.view.remove_component(self.current_component)
        mol = self.df.iloc[self.frm.value]['ROMol']
        mol_comb = self.merge_mols(lig_mol=mol)
        fh = StringIO(Chem.MolToPDBBlock(mol_comb))
        self.current_component = self.view.add_component(fh, ext='pdb')
        self.current_component.center(selection='hetero', duration=1000)
        self.current_component.add_representation('contact', selection='*', weak_hydrogen_bond=True)
        self.view.add_representation('line', selection='protein')

        return self.view

    def box(self, chg=None):
        self.frm = BoundedIntText(
            value=0,
            min=0,
            max= len(self.df),
            step=1,
            description='Display mol:',
            disabled=False
        )

        return self.frm

    def select_button(self, chg=None):
        self.select_button = Button(
            description='Display',
            disabled=False,
            button_style='success',
            tooltip='Display',
            icon='check'
        )
        self.select_button.on_click(self.add_mol)

        return self.select_button

    def show_prot_surface(self, chg=None):

        self.prot_surface = self.current_component.add_representation('surface',
                                                                      selection='protein',
                                                                      color_scheme='electrostatic',
                                                                      surface_type='av',
                                                                      opacity=0.6)

    def show_lig_surface(self, chg=None):

        self.prot_surface = self.current_component.add_representation('surface',
                                                                      selection='hetero',
                                                                      color_scheme='electrostatic',
                                                                      surface_type='av',
                                                                      opacity=0.6)

    def surface_button(self, chg=None):
        self.surface_button = Button(
            description='Show prot surface',
            disabled=False,
            button_style='info',
            tooltip='Show surface',
        )
        self.surface_button.on_click(self.show_prot_surface)

        return self.surface_button

    def lig_surface_button(self, chg=None):
        self.surface_button = Button(
            description='Show lig surface',
            disabled=False,
            button_style='warning',
            tooltip='Show surface',
        )
        self.surface_button.on_click(self.show_lig_surface)

        return self.surface_button

    def viewer(self):
        return self.view
