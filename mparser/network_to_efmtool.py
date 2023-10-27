# -*- coding: utf-8 -*-

"""
network_to_efmtool.py

Converts a metabolism network to EFMTool files

"""

from .meta_network import MetaNetwork
from .parsehelper import add_quotes

def to_efmtool(reactions, reversibles, metabolites, metint, stoichiometry, matrix):
    """
    Converts every list created from reading the catalyzed reactions to EFMTOOL input format
    Excludes external metabolites and their stoichiometry
    
    Params:
        reactions: reactions set
        reversibles: reversible reactions set
        metabolites: metabolites set
        metint: set of internal metabolites
        stoichiometry: list of stoichiometry tuples (metabolite name, reaction name, value)
        matrix: stoichiometry matrix
        
    Returns:
        fcontent: list of reactions with reaction names, metabolite names and stoichiometry
                  file formats accepted as input by EFMTool             
    """
    efmtool = {'mnames.txt': "",
               'revs.txt': "",
               'rlist.txt': "",
               'rnames.txt': "",
               'stoich.txt': ""} # NotImplemented
    
    for metabolite in metint:
        efmtool['mnames.txt'] += add_quotes(metabolite) + '\t'
        
    for reaction in reactions:
        reversible = reaction in reversibles
        efmtool['rnames.txt'] += add_quotes(reaction) + '\t'
        efmtool['revs.txt'] += str(int(reversible)) + '\t'
        metatxt = lambda meta, coef: str(abs(coef)) + ' ' + meta if abs(coef) != 1 else meta 
        reacts = [metatxt(x, z) for x, y, z in stoichiometry if y == reaction and z < 0 and x in metint] # uncool
        prods = [metatxt(x, z) for x, y, z in stoichiometry if y == reaction and z > 0 and x in metint]
        arrow = '<-->' if reversible else '-->'
        formula = ' + '.join(reacts) + ' ' + arrow + ' ' + ' + '.join(prods)
        efmtool['rlist.txt'] += add_quotes(reaction) + '\t'+ add_quotes(formula) + '\n'
        
    efmtool['stoich.txt'] += str(matrix) 
        
    fcontent = ""
    for fname in efmtool:
        fcontent += f"---- {fname} ----\n\n"
        fcontent += efmtool[fname].strip() + '\n\n'
        
    return fcontent      
    

def format_efmtool(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to EFMTool input format
    
    Params:
        network: MetaNetwork object
        out_fname: EFMTool output file
    """
    efmtool = to_efmtool(network.reactions, network.reversibles, network.metabolites,
                         network.metint, network.stoichiometry, network.matrix)
    f_out = open(out_fname, "w")
    f_out.write(efmtool)
    f_out.close()  
    
