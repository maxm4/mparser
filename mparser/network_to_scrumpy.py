# -*- coding: utf-8 -*-

"""
network_to_scrumpy.py

Converts a metabolism network to ScrumPy SPY format

"""

from .meta_network import MetaNetwork
from .parsehelper import add_quotes
import numpy as np

def to_scrumpy(reactions, reversibles, metabolites, metint, metext, stoichiometry):
    """
    Converts every list created from reading the catalyzed reactions to ScrumPy input format    
    Params:
        reactions: reactions set
        reversibles: reversible reactions set
        metabolites: metabolites set
        metint: set of internal metabolites
        metext: set of external metabolites
        stoichiometry: list of stoichiometry tuples (metabolite name, reaction name, value)
        
    Returns:
        fcontent: text file format accepted by ScrumPy             
    """
    
    fcontent="Structural()\n"
        
    for reaction in reactions:
        fcontent += '\n'
        reversible = reaction in reversibles
        disp_coef = lambda coef: np.format_float_positional(abs(coef), trim='-')
        auto_x = lambda meta: 'x_' + meta if meta in metext else meta
        metatxt = lambda meta, coef: disp_coef(coef) + ' ' + auto_x(meta) if abs(coef) != 1 else auto_x(meta) 
        reacts = [metatxt(x, z) for x, y, z in stoichiometry if y == reaction and z < 0]
        prods = [metatxt(x, z) for x, y, z in stoichiometry if y == reaction and z > 0]
        arrow = '<>' if reversible else '->'
        formula = ' + '.join(reacts) + ' ' + arrow + ' ' + ' + '.join(prods)
        fcontent += reaction + ':\n'
        fcontent += formula + ' ~\n'
        
    return fcontent      
    

def format_scrumpy(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to ScrumPy input format
    
    Params:
        network: MetaNetwork object
        out_fname: ScrumPy output file
    """
    efmtool = to_scrumpy(network.reactions, network.reversibles, network.metabolites,
                         network.metint, network.metext, network.stoichiometry)
    f_out = open(out_fname, "w")
    f_out.write(efmtool)
    f_out.close()  
