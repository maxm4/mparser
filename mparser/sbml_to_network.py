# -*- coding: utf-8 -*-

"""
sbml_to_network.py

Reads a metabolism network in SBML format
Returns a MetaNetwork object

"""

from . import parsehelper as ph
from .meta_network import MetaNetwork
from libsbml import readSBMLFromString
from .global_parameters import *

def is_null_reaction(reaction, metext):
    """
    Checks if a reaction consumes and produces at least one internal metabolite
    
    Params:
        reaction: libSBML model Reaction object
        metext: list of external metabolites
        
    Returns:
        boolean, true if the reaction does not consume or produce any internal metabolite
    """
    if reaction.getNumReactants() == 0 or reaction.getNumProducts() == 0:
        return True
    for j in range(reaction.getNumReactants()):
        if reaction.getReactant(j).getSpecies() not in metext:
            return False
    for k in range(reaction.getNumProducts()):
        if reaction.getProduct(k).getSpecies() not in metext:
            return False
    return True



def read_sbml(f_in, metext_func=metext_func):
    """
    Reads network from SBML file
    Converts it to a MetaNetwork object
    
    Params:
        f_in: input metabolism network file descriptor
        metext_func: function determining whether if a metabolite is external
            may be one of the lambda functions defined in this module:
                boundary_condition, external_compartment, etc.
        
    Returns:
        network: MetaNetwork object
    """
    strfile = f_in.read()
    document = readSBMLFromString(strfile)
    
    if document.getNumErrors() > 0:
        document.printErrors()
        raise ValueError("Encountered the SBML errors listed above")

    model = document.getModel()

    if model is None:
        raise ValueError("No model present")

    species = [model.getSpecies(i) for i in range(model.getNumSpecies())]
    metabolites = [specie.getId() for specie in species]
    metext = [specie.getId() for specie in species if metext_func(specie)]    

    reactions = []
    reversibles = []
    stoichiometry = []
    nullreacs = []

    for i in range(0,model.getNumReactions()):
        r = model.getReaction(i)
        reac = r.getId()
            
        if is_null_reaction(r, metext):
            if remove_null_reactions:
                nullreacs.append(reac)
                continue # exclude reaction

        reactions.append(reac)
        if r.getReversible():
            reversibles.append(reac)

        for j in range(0,r.getNumReactants()):
            rp = r.getReactant(j)
            coeff = -rp.getStoichiometry()
            meta = rp.getSpecies()
            stoichiometry.append((meta, reac, coeff))

        for j in range(0,r.getNumProducts()):
            rp = r.getProduct(j)
            coeff = rp.getStoichiometry()
            meta = rp.getSpecies()
            stoichiometry.append((meta, reac, coeff))

    return MetaNetwork(metext, metabolites, reactions, stoichiometry, reversibles, nullreacs)    

def sbml_read_network(in_fname):
    """
    Reads a .txt metabolism network in MetaTool format and returns a MetaNetwork object
    
    Params:
        in_fname: input metabolism network file name

    Returns:
        network: MetaNetwork object
    """
    f_in = open(in_fname)
    network = read_sbml(f_in)
    f_in.close()
    return network
