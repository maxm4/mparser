# -*- coding: utf-8 -*-

"""
network_to_sbml.py

Converts a metabolism network to a SBML file

"""

from .meta_network import MetaNetwork
from libsbml import SBMLNamespaces, SBMLDocument, writeSBMLToString, writeSBMLToFile
import re, os

def to_sbml(reactions, reversibles, metabolites, metext, metint, stoichiometry, name, fbc_plugin=True):
    """
    Converts every list created from reading the catalyzed reactions to a SBML file
    
    Params:
        reactions: reactions list
        reversibles: reversible reaction set
        metabolites: metabolites list
        metext: set of external metabolites
        metint: set of internal metabolites
        stoichiometry: list of stoichiometry tuples (metabolite name, reaction name, value)
        name: name of the model
        fbc_plugin: boolean indicating whether the SBML file should include the fbc plugin
        
    Returns:
        sbml: sbml format string
    """
    
    try:
        # SBML Level 3 Version 2 with fbc plugin Version 2
        sbmlns = SBMLNamespaces(3,2,"fbc", 2) if fbc_plugin else SBMLNamespaces(3,2)
        document = SBMLDocument(sbmlns)
    except ValueError:
        raise SystemExit('Could not create SBMLDocument object')
        
    model = document.createModel()

    if fbc_plugin:
        document.setPackageRequired("fbc", False)  
        mplugin = model.getPlugin("fbc")
        
    # SBML ids must be alphanumeric, '_' is the only special character allowed
    # Matches special characters and replaces them with '_'
    sbml_id = lambda s: re.sub(r'\W', '_', s)

    # Create an external compartment with external metabolites

    cext = model.createCompartment()
    cext.setId('ext') 
    cext.setConstant(True)
    
    for met in metext:
        s1 = model.createSpecies()
        s1.setName(met)
        s1.setId('M_' + sbml_id(met))
        s1.setCompartment('ext')
        s1.setBoundaryCondition(True) # does not appear in the matrix
        s1.setConstant(False)
        s1.setHasOnlySubstanceUnits(False)
        
    # Create an internal compartment with internal metabolites

    cint = model.createCompartment()
    cint.setId('int')
    cint.setConstant(True)
    
    for metabolite in metint:
        s1 = model.createSpecies()
        s1.setName(metabolite)
        s1.setId('M_' + sbml_id(metabolite))
        s1.setCompartment('int')
        s1.setBoundaryCondition(False) # does appear in the matrix
        s1.setConstant(False)
        s1.setHasOnlySubstanceUnits(False)
        
    
    # Create flux bounds 
    
    if fbc_plugin:
        bound = 1000
    
        lb = model.createParameter()
        lb.setId('lb')
        lb.setName("lower bound")
        lb.setConstant(True)
        lb.setValue(-bound)
    
        lb = model.createParameter()
        lb.setId('zb')
        lb.setName("zero bound")
        lb.setConstant(True)
        lb.setValue(0)
        
        lb = model.createParameter()
        lb.setId('ub')
        lb.setName("upper bound")
        lb.setConstant(True)
        lb.setValue(bound)
            
    # Add reactions

    for reaction in reactions:
        reversible = reaction in reversibles
        r1 = model.createReaction()
        r1.setName(reaction)
        r1.setId('R_' + sbml_id(reaction))
        r1.setFast(False)
        r1.setReversible(reversible)
        if fbc_plugin:
            f1 = r1.getPlugin('fbc')
            f1.setLowerFluxBound('lb' if reversible else 'zb')
            f1.setUpperFluxBound('ub')
        
    # Add reactants and products - stoichiometry  
    
    createSub = (lambda r, x: r.createReactant if x < 0 else r.createProduct)
    for stuple in stoichiometry:
        reaction = model.getReaction('R_' + sbml_id(stuple[1]))
        rp = createSub(reaction, stuple[2])()
        rp.setSpecies('M_' + sbml_id(stuple[0]))
        rp.setStoichiometry(abs(stuple[2]))
        rp.setConstant(True)

    model.setName(name)
    model.setId(sbml_id(name))
    model.setMetaId(sbml_id(name))

    # Set objective
    if fbc_plugin:
        objname = reactions[0]
        biomass = [s for s in reactions if 'biomass' in s.lower()]
        if biomass:
            objname = biomass[0]  
        objective = mplugin.createObjective()
        objective.setId("obj")
        objective.setType("maximize")
        maxflux = objective.createFluxObjective()
        maxflux.setReaction('R_' + sbml_id(objname))
        maxflux.setCoefficient(1)
        mplugin.setActiveObjectiveId("obj")

    mplugin.setStrict(True)

    # Return a text string containing the model in XML format.

    return writeSBMLToString(document)       
    

def format_sbml(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to SBML
    
    Params:
        network: MetaNetwork object
        out_fname: output SBML file name
    """
    sbml = to_sbml(network.reactions, network.reversibles, network.metabolites,
                   sorted(network.metext), network.metint, network.stoichiometry, os.path.basename(out_fname))
    f_out = open(out_fname, "w")
    f_out.write(sbml)
    f_out.close()  
    
