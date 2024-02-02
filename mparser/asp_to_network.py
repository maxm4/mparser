# -*- coding: utf-8 -*-

"""
asp_to_network.py

Reads a metabolism network - in ASP format
Returns a MetaNetwork object

"""

from . import parsehelper as ph
from .meta_network import MetaNetwork 

get_inside_pars = lambda s: s[s.find("(")+1:s.find(")")]
   
def ground_asp(f_in):
    """
    Grounds ASP atoms
    
    Params:
        f_in: input metabolism network file descriptor
        
    Returns:
        lines: lines read
    """      
    lines = []
    for s in f_in.readlines():
        deb, inside, end = s[:s.find("(")+1], s[s.find("(")+1:s.find(")")], s[s.find(")"):]
        insides = inside.split(';')
        for content in insides:
            line = deb + content + end
            lines.append(line)
    #print(lines)
    return lines
            
       
def read_asp(lines):
    """
    Reads ASP atoms
    Converts it to a MetaNetwork object
    
    Params:
        lines: Lines of the ASP files
        
    Returns:
        network: MetaNetwork object
    """    
    reactions = []
    reversibles = []
    externals = []
    metabolites = []
    stoichiometry = []
    for line in lines:
        if line.startswith('reaction'):
            reaction = eval(get_inside_pars(line))
            if reaction.endswith('_rev'):
                reversibles.append(reaction[:-4])
            else:
                reactions.append(reaction)
        elif line.startswith('external'):
            external = eval(get_inside_pars(line))
            externals.append(external)
        elif line.startswith('metabolite'):
            metabolite = eval(get_inside_pars(line))
            metabolites.append(metabolite)
        elif line.startswith('stoichiometry'):
            stoich = eval(get_inside_pars(line))                   
            assert(len(stoich) == 3)            
            if not stoich[1].endswith('_rev'):
                stoichiometry.append((str(stoich[0]), str(stoich[1]), float(stoich[2])))
    metabolites.extend(externals)
    externals = set(externals)
    metabolites = set(metabolites)
    reactions = set(reactions)
    stoichiometry = set(stoichiometry)
    reversibles = set(reversibles)
    return MetaNetwork(externals, metabolites, reactions, stoichiometry, reversibles)
    
    

def asp_read_network(in_fname):
    """
    Reads a .txt metabolism network in MetaTool format and returns a MetaNetwork object
    
    Params:
        in_fname: input metabolism network file name

    Returns:
        network: MetaNetwork object
    """
    f_in = open(in_fname)
    lines = ground_asp(f_in)
    network = read_asp(lines) 
    f_in.close()
    return network
