import sys
import os
import ipywidgets as ui
from IPython.display import display


class PathSelector():

    def __init__(self, start_dir, select_file=True):
        self.file = None
        self.select_file = select_file
        self.cwd = start_dir
        self.select = ui.SelectMultiple(options=['init'], value=(), rows=10, description='')
        self.accord = ui.Accordion(children=[self.select])

        self.accord.selected_index = None # Start closed (showing path only)
        self.refresh(self.cwd)
        self.select.observe(self.on_update,'value')

    def on_update(self, change):
        if len(change['new']) > 0:
            self.refresh(change['new'][0])

    def refresh(self, item):
        path = os.path.abspath(os.path.join(self.cwd, item))

        if os.path.isfile(path):
            if self.select_file:
                self.accord.set_title(0, path)
                self.file = path
                self.accord.selected_index = None
            else:
                self.select.value = ()

        else: # os.path.isdir(path)
            self.file = None
            self.cwd  = path

            # Build list of files and dirs
            keys = ['[..]']
            for item in os.listdir(path):
                if item[0] == '.':
                    continue
                elif os.path.isdir(os.path.join(path,item)):
                    keys.append('['+item+']')
                else:
                    keys.append(item)

            # Sort and create list of output values
            keys.sort(key=str.lower)
            vals = []
            for k in keys:
                if k[0] == '[':
                    vals.append(k[1:-1]) # strip off brackets
                else:
                    vals.append(k)

            # Update widget
            self.accord.set_title(0,path)
            self.select.options = list(zip(keys, vals))
            with self.select.hold_trait_notifications():
                self.select.value = ()

class DockButtons():
    def __init__(self, path_selector):
        
        item_layout = ui.Layout(margin='0px 50px 0px 5px')
        
        self.f = path_selector
        self.pdb_button = ui.Button(description="Set PDB file", layout=item_layout)
        self.pdb_button.on_click(self.on_button_clicked)

        self.lig_button = ui.Button(description="Set ligand(s) SDF file", layout=item_layout)
        self.lig_button.on_click(self.on_button_clicked)

        self.ref_button = ui.Button(description="Set ref ligand SDF file", layout=item_layout)
        self.ref_button.on_click(self.on_button_clicked)
        
        
    
    def fun(self, a):
        a = self.f.file
        return a

    def on_button_clicked(self, b):
        pth = self.fun(b)
        if not pth:
            print('PLEASE SELECT A PATH FROM THE THING ABOVE!')
            return
        b.description = pth
        
    def render_buttons(self):

        widge = ui.HBox([ui.Label('Protein:'), self.pdb_button, 
                      ui.Label('Query Ligands:'), self.lig_button, 
                      ui.Label('Reference Ligand:'), self.ref_button])

    def render_one_button(self, label):
        button = ui.Button(description=label)
        button.on_click(self.on_button_clicked)

    
        return button


class InputBoxSet():
    def __init__(self, tlabel):
        self.inbox = ui.Text(description=tlabel)
        self.button = ui.Button(description='Set', icon='check')
        self.button.on_click(self.on_button_clicked)

        self.button2 = ui.Button(description='Create', icon='plus-square')
        self.button2.on_click(self.create_dir)

    def fun(self, a):
        if os.path.isdir(self.inbox.value):
            a = self.inbox.value
            print(a)
            return a
        else:
            raise Exception('please input a valid path!')

    def on_button_clicked(self, b):
        pth = self.fun(b)
        b.description = pth

    def create_dir(self, b):
        if not os.path.isdir(self.inbox.value):
            os.makedirs(self.inbox.value)
            print(f'Created: {self.inbox.value}')
        else:
            print(f'Directory {self.inbox.value} exists!')

    def render(self):

        widge = ui.HBox([self.inbox, self.button, self.button2])

        return widge

