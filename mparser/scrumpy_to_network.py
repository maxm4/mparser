# -*- coding: utf-8 -*-

"""
scrumpy_to_network.py

Reads a metabolism network
Returns a MetaNetwork object

"""

from . import parsehelper as ph
from .meta_network import MetaNetwork 

qtreat = ph.smart_remove_quotes 
    
def read_metext(f_in):
    """
    Reads external metabolites from ScrumPy file
    
    Params:
        f_in: input metabolism network file descriptor
        
    Returns:
        metext: list of external metabolites
    """
    line = f_in.readline() # skip comments and whitespace
    while line.startswith('#') or line.isspace():
        line = f_in.readline()
    line = line.strip()
    assert(line == "Structural()") # file must begin by "Structural"
    line = f_in.readline().strip()
    if line.startswith("External"):
       return line[9:-1].split(', ')
    return []
    
    
def read_cat(f_in, metext):
    """
    Reads catalyzers from ScrumPy file
    Converts it to a MetaNetwork object
    
    Params:
        f_in: input metabolism network file descriptor
        metext: list of external metabolites
        
    Returns:
        network: MetaNetwork object
    """
    metext = set(metext)
    reactions = []
    reversibles = []
    metaboset = set()
    metabolites = []
    stoichiometry = []
    for line in f_in.readlines():
        if line.isspace():
            continue
        if line.startswith("#"): # eliminates empty lines and comments
            continue
        if not ':' in line:
            line = line.strip()
            line = line.replace('->', '-> ')
            line = line.replace('<>', '<> ')
            line = line.replace('<-', '<- ') 
            if line.startswith('~'):
                continue
            if line.endswith('~'):
                line = line[:-1]
            words = line.split()
            pos_i = ph.indexl(words, '->') # prioritize this list index if not -1
            pos_r = ph.indexl(words, '<>')
            pos_b = ph.indexl(words, '<-')
            #print("###############", words)
            assert(len(words) > 2)  
            assert(pos_i != -1 or pos_r != -1 or pos_b != -1)
            
            backwards = 1
            if pos_i == -1 and pos_b == -1:
                reversibles.append(reaction)
                pos = pos_r
            elif pos_i == -1:
                backwards = -1
                pos = pos_b
            else:
                pos = pos_i
                
            coeff = 1
            for j in range(0, len(words)):
                if words[j] == "+" or j == pos:
                    coeff = 1
                elif ph.iswint(words[j]):
                    coeff = int(words[j])
                elif ph.iswfloat(words[j]):
                    coeff = float(words[j])
                else:
                    metabolite = qtreat(words[j])
                    if metabolite.startswith("x_"):
                        metext.add(metabolite)
                    if metabolite not in metaboset:
                        metaboset.add(metabolite)
                        metabolites.append(metabolite)
                    if j < pos:
                        coeff = -coeff
                    coeff *= backwards
                    stoichiometry.append((metabolite, reaction, coeff))
        else:
            line = line.strip()
            words = line.split(':')
            reaction = qtreat(words[0])
            reactions.append(reaction)

    
    return MetaNetwork(metext, metabolites, reactions, stoichiometry, reversibles)
    
    

def scrumpy_read_network(in_fname):
    """
    Reads a .spy metabolism network in ScrumPy format and returns a MetaNetwork object
    
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
