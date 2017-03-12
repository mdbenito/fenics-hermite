# Imports IPython notebooks.
# Based on https://github.com/adrn/ipython/blob/master/examples/Notebook/Importing%20Notebooks.ipynb
# See also: https://github.com/grst/nbimporter
#
# This basically adds some parsing which avoids executing any code
# in the notebooks which isn't a definition. See NotebookLoader below.

import io, os, sys, types, ast
import nbformat

from IPython import get_ipython
from IPython.core.interactiveshell import InteractiveShell

def find_notebook(fullname, path=None):
    """ Find a notebook, given its fully qualified name and an
    optional path.
    
    This turns "foo.bar" into "foo/bar.ipynb" and tries turning 
    "Foo_Bar" into "Foo Bar" if Foo_Bar does not exist.
    """
    name = fullname.rsplit('.', 1)[-1]
    if not path:
        path = ['']
    for d in path:
        nb_path = os.path.join(d, name + ".ipynb")
        if os.path.isfile(nb_path):
            return nb_path
        # let import Notebook_Name find "Notebook Name.ipynb"
        nb_path = nb_path.replace("_", " ")
        if os.path.isfile(nb_path):
            return nb_path

class CellTransformer(ast.NodeTransformer):
    """ Removes all nodes from an AST tree which are not suitable
    for exporting out of a notebook. """
    def visit(self, node):
        """ Visit a node. """
        if node.__class__.__name__ in ['Module', 'FunctionDef', 'ClassDef', 
                                       'Import', 'ImportFrom']:
            return node
        return None

class NotebookLoader(object):
    """ Module Loader for IPython Notebooks.
    
    Importing from a notebook is different from a module: because one
    typically has many computations and tests besides actual usable code,
    here we only run code which either defines a function or a class, or
    imports code from other modules and notebooks.
    
    Furthermore, in order to provide per-notebook initialisation, if a
    special function __init__() is defined in the notebook, it will be
    executed the first time an import statement is.
    """
    
    def __init__(self, path=None):
        self.shell = InteractiveShell.instance()
        self.path = path
    
    def load_module(self, fullname):
        """ Imports a notebook as a module. """
        path = find_notebook(fullname, self.path)
        
        #print("Importing notebook `%s`." % path)
                                       
        # load the notebook object
        nb = nbformat.read(path, as_version=4)    
        
        # create the module and add it to sys.modules
        # if name in sys.modules:
        #    return sys.modules[name]
        mod = types.ModuleType(fullname)
        mod.__file__ = path
        mod.__loader__ = self
        mod.__dict__['get_ipython'] = get_ipython
        sys.modules[fullname] = mod
        
        # extra work to ensure that magics that would affect the user_ns
        # actually affect the notebook module's ns
        save_user_ns = self.shell.user_ns
        self.shell.user_ns = mod.__dict__
        
        try:
            deleter = CellTransformer()
            for cell in filter(lambda c: c.cell_type == 'code', nb.cells):
                # transform the input to executable Python
                code = self.shell.input_transformer_manager.transform_cell(cell.source)
                # Remove anything that isn't a def or a class
                tree = deleter.generic_visit(ast.parse(code))
                # run the code in the module
                codeobj = compile(tree, filename=path, mode='exec')
                exec(codeobj, mod.__dict__)
        finally:
            self.shell.user_ns = save_user_ns

        # Run any initialisation if available, but only once
        if not mod.__dict__.has_key('__init_done__'):
            try:
                mod.__dict__['__init__']()
                mod.__init_done__ = True
            except KeyError:
                pass

        return mod

class NotebookFinder(object):
    """ Module finder that locates IPython Notebooks. """
    def __init__(self):
        self.loaders = {}
    
    def find_module(self, fullname, path=None):
        nb_path = find_notebook(fullname, path)
        if not nb_path:
            return
        
        key = path
        if path:
            # lists aren't hashable
            key = os.path.sep.join(path)
        
        if key not in self.loaders:
            self.loaders[key] = NotebookLoader(path)
        return self.loaders[key]


# Import stuff
sys.meta_path.append(NotebookFinder())
