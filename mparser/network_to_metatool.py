# -*- coding: utf-8 -*-

"""
network_to_metatool.py

Converts a metabolism network to MetaTool TXT format

"""

from .meta_network import MetaNetwork
from .parsehelper import add_quotes
import warnings

def to_metatool(reactions, reversibles, metabolites, metint, metext, stoichiometry):
    """
    Converts every list created from reading the catalyzed reactions to MetaTool input format    
    Params:
        reactions: reactions set
        reversibles: reversible reactions set
        metabolites: metabolites set
        metint: set of internal metabolites
        metext: set of external metabolites
        stoichiometry: list of stoichiometry tuples (metabolite name, reaction name, value)
        
    Returns:
        fcontent: text file format accepted by MetaTool             
    """
    
    fcontent="-METEXT\n"
    for met in metext:
        fcontent += met + ' '
    fcontent = fcontent[:-1] + '\n'
    fcontent += '\n'
    
    fcontent += '-CAT\n'     
    for reaction in reactions:
        reversible = reaction in reversibles
        metatxt = lambda meta, coef: str(abs(coef)) + ' ' + meta if abs(coef) != 1 else meta 
        reacts = [metatxt(x, z) for x, y, z in stoichiometry if y == reaction and z < 0]
        prods = [metatxt(x, z) for x, y, z in stoichiometry if y == reaction and z > 0]
        arrow = '=' if reversible else '=>'
        formula = ' + '.join(reacts) + ' ' + arrow + ' ' + ' + '.join(prods)
        if (not reacts) or (not prods):
            warnings.warn('Producing a file that is unreadable by MetaTool due to reactions with only production/consumption') 
            warnings.warn('If this is the result of conversion of a Cobra model, please run conversion on a model with boundaries created by cobra_to_sbml.py instead') 
        fcontent += reaction + ' : '+ formula + '\n'
        
    return fcontent      
    

def format_metatool(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to MetaTool input format
    
    Params:
        network: MetaNetwork object
        out_fname: MetaTool output file
    """
    efmtool = to_metatool(network.reactions, network.reversibles, network.metabolites,
                         network.metint, network.metext, network.stoichiometry)
    f_out = open(out_fname, "w")
    f_out.write(efmtool)
    f_out.close()  
