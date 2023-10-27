# -*- coding: utf-8 -*-

"""
metatool_to_network.py

Reads a metabolism network
Returns a MetaNetwork object

"""

from . import parsehelper as ph
from .meta_network import MetaNetwork 
    
def read_metext(f_in):
    """
    Reads external metabolites from MetaTool .txt file
    
    Params:
        f_in: input metabolism network file descriptor
        
    Returns:
        metext: list of external metabolites
    """
    line = f_in.readline().strip()
    assert(line == "-METEXT") # file must begin by "-METEXT"
    line = f_in.readline().strip()      
    return line.split() 
    
    
def read_cat(f_in, metext):
    """
    Reads catalyzers from MetaTool .txt file
    Converts it to a MetaNetwork object
    
    Params:
        f_in: input metabolism network file descriptor
        metext: list of external metabolites
        
    Returns:
        network: MetaNetwork object
    """
    
    line = f_in.readline().strip()
    if line != "-CAT":
        line = f_in.readline().strip()
    assert(line == "-CAT") # that line of the file must be "-CAT"
    
    reactions = []
    reversibles = []
    metaboset = set()
    metabolites = []
    stoichiometry = []
    for line in f_in.readlines():
        line = line.strip()
        if line and not line.startswith("#"): # eliminates empty lines and comments
            words = line.split()
            pos_i = ph.indexl(words, '=>') # prioritize this list index if not -1
            pos_r = ph.indexl(words, '=')
            
            assert(len(words) > 4)
            assert(words[1] == ":")   
            assert(pos_i != -1 or pos_r != -1)
            
            reaction = words[0] 
            reactions.append(reaction)
            
            if pos_i == -1:
                reversibles.append(reaction)
                pos = pos_r
            else:
                pos = pos_i
                
            coeff = 1
            for j in range(2, len(words)):
                if words[j] == "+" or j == pos:
                    coeff = 1
                elif ph.iswint(words[j]):
                    coeff = int(words[j])
                elif ph.iswfloat(words[j]):
                    coeff = float(words[j])
                else:
                    metabolite = words[j]
                    if metabolite not in metaboset:
                        metaboset.add(metabolite)
                        metabolites.append(metabolite)
                    if j < pos:
                        coeff = -coeff
                    stoichiometry.append((metabolite, reaction, coeff))
    
    return MetaNetwork(metext, metabolites, reactions, stoichiometry, reversibles)
    
    

def txt_read_network(in_fname):
    """
    Reads a .txt metabolism network in MetaTool format and returns a MetaNetwork object
    
    Params:
        in_fname: input metabolism network file name

    Returns:
        network: MetaNetwork object
    """
    f_in = open(in_fname)
    metext = read_metext(f_in)
    network = read_cat(f_in, metext) 
    f_in.close()
    return network
