# -*- coding: utf-8 -*-

"""
network_to_asp.py

Converts a metabolism network to ASP rules

"""

from .meta_network import MetaNetwork
from .parsehelper import add_quotes as aq

capitalize_plural = lambda x: x[0].upper() + x[1:] + 's'

def to_asp(reactions, reversibles, metabolites, metext, extras, stoichiometry, reverse_stoichiometry):
    """
    Converts every list created from parsing the network to ASP rules
    
    Params:
        reactions: reactions list
        reversibles: reversible reactions set
        metabolites: metabolites list
        metext: set of external metabolites
        extras: extra fields dict of lists
        stoichiometry: string containing stoichiometry ASP constraints
        reverse_stoichiometry:  string containing stoichiometry ASP constraints
                                from reversible reactions
        
    Returns:
        asp: ASP code string
    """
    
    asp = str()
    if metext:
        asp += "% External metabolites \n"
        asp += "external("
        for met in metext: # would have been nice to have the list here to keep the reading order
            asp += aq(met) + ";"
        asp = asp[:-1] + ").\n"  
        
    asp += "% Internal metabolites \n"
    asp += "metabolite("
    for metabolite in metabolites:
        if metabolite not in metext:
            asp += aq(metabolite) + ";"
    asp = asp[:-1] + ").\n" 
    
    if reversibles:    
        asp += "% Reversible reactions \n"
        asp += "reversible("
        for reaction in reversibles:
            reaction_rev = reaction + "_rev"
            asp += aq(reaction) + "," + aq(reaction_rev) + ";"
        asp = asp[:-1] + ").\n"  
        
    asp += "% All reactions \n"
    asp += "reaction("
    for reaction in reactions:
        asp += aq(reaction) + ";"
        if reaction in reversibles:
            reaction_rev = reaction + "_rev"
            asp += aq(reaction_rev) + ";"
    asp = asp[:-1] + ").\n"

    for extra, field in extras.items():
        if field:
            asp += f"% {capitalize_plural(extra)} \n"
            asp += f"{extra}("
            for element in field:
                asp += aq(element) + ";"
                if (element in reactions) and (element in reversibles):
                    reaction_rev = element + "_rev"
                    asp += aq(reaction_rev) + ";"
            asp = asp[:-1] + ").\n" 
            
    asp += "% Stoichiometry \n"
    asp += stoichiometry  
    
    asp += "\n% Reversible reactions stoichiometry \n"
    asp += reverse_stoichiometry  

    return asp         
    
    
    
def make_stoichiometry(metabo, reaction, coeff, inv=False):
    """
    Makes an ASP string describing the stoichiometry
    
    Params:
        metabo: metabolite name
        reaction: reaction name
        coeff: stoichiometric coefficient
        inv: if True, computes inverse of the reaction
        
    Returns:
        asp: ASP code string
    """    
    if not isinstance(coeff, float):
        coeff = float(coeff)

    if inv:
        reaction += "_rev"
        coeff = -coeff

    if isinstance(coeff, float):
        coeff = aq(str(coeff))

    asp = "stoichiometry(" + aq(metabo) + "," + aq(reaction) + "," + str(coeff) + ").\n"
    return asp    


def make_asp(network):
    """
    Makes ASP constraints for the metabolism network
    
    Params:
        network: MetaNetwork object
    """
    stoichiometry = ""
    reverse_stoichiometry = ""
    for stuple in network.stoichiometry:
        stoichiometry += make_stoichiometry(stuple[0], stuple[1], stuple[2])
        if stuple[1] in network.reversibles:
            reverse_stoichiometry += make_stoichiometry(stuple[0], stuple[1], stuple[2], inv=True)
    return to_asp(network.reactions, network.reversibles, network.metabolites, network.metext,
                  network.extras, stoichiometry, reverse_stoichiometry)
    

def format_asp(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to ASP
    
    Params:
        network: MetaNetwork object
        out_fname: output ASP file name
    """
    asp = make_asp(network)
    f_out = open(out_fname, "w")
    f_out.write(asp)
    f_out.close()  
    
    
def write_asp(asp:str, out_fname):
    """
    Writes an ASP string to a file
    
    Params:
        asp: ASP string
        out_fname: output ASP file name
    """
    f_out = open(out_fname, "w")
    f_out.write(asp)
    f_out.close()  
    

def append_asp(asp:str, out_fname):
    """
    Appends an ASP string to a file
    
    Params:
        asp: ASP string
        out_fname: output ASP file name
    """
    f_out = open(out_fname, "a")
    f_out.write(asp)
    f_out.close() 
    
