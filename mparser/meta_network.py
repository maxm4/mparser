# -*- coding: utf-8 -*-

"""
meta_network.py

Module containing the MetaNetwork class
Used to store metabolic network information when reading and converting to ASP

"""

import warnings
import numpy as np

class MetaNetwork:
    """
    
    MetaNetwork class
    Class for storing metabolic network information
    
    """ 

    def only_one_stoch_check(self):
        """
        Checks for reactions using a metabolite as both a substrate and product
        A given metabolite and reaction should be associated to only one stoichiometry coefficient
        
        """
        only_one = {}
        for metabolite, reaction, _ in self.stoichiometry:
            if only_one.get((metabolite, reaction)) is not None:
               raise ValueError(f'Metabolite {metabolite} should be only used once in reaction {reaction}') 
            only_one[(metabolite, reaction)] = True
        
    
    def null_reactions_check(self):
        """
        Checks for reactions with no internal metabolite production or consumption
            
        """
        for reaction in self.reactions:
            metas = [stoch[0] for stoch in self.stoichiometry if stoch[1] == reaction]
            if len(metas) < 2:
                warnings.warn('Found reactions with missing overall metabolite production or consumption')
            elif set(metas).issubset(self.metext):
                warnings.warn('Found reactions with missing internal metabolite production or consumption')
        
    
    def structures_check(self):
        """
        Checks the consistency of the network data structures.
        Raises an error if an incoherence is found.
            
        """
        if not self.metabolites:
            raise ValueError('No metabolite')
        if not self.reactions:
            raise ValueError('No reaction')
        if not len(self.metabolites) == len(self._metabolites_s):
            raise ValueError('Duplicate found in metabolites list')
        if not len(self.reactions) == len(self._reactions_s):
            raise ValueError('Duplicate found in reactions list')
        if not len(self.reversibles) == len(self._reversibles_l):
            raise ValueError('Duplicate found in reversible reactions list')
        if not len(self.metext) == len(self._metext_l):
            raise ValueError('Duplicate found in external metabolites list')
        if not self.metext:
            warnings.warn('No external metabolite definition', stacklevel=2)
        if self.nb_nullreacs > 0:
            warnings.warn(f"{self.nb_nullreacs} incorrect reactions were excluded from the" +
                          " model with the current external metabolite definition", stacklevel=2)
        if not self.metext.issubset(self._metabolites_s):
            print(self.metext - self._metabolites_s) 
            raise ValueError('At least one external metabolite was not found in the metabolite list')
        if not self.reversibles.issubset(self._reactions_s):
            raise ValueError('At least one reversible reaction was not found in the reaction list')
            
            
    def compute_stats(self):
        """
        Computes some caracteristics of the model, such as the number of external metabolites.
           
        Initializes:
            nb_metext: number of external metabolites
            nb_transporters: number of transport or exchange reactions
            nb_reversibles: number of reversible reactions
            nb_nullreacs: number of exchange reactions excluded from the model
        
        """
        self.nb_metext = len(self._metext_l)
        self.nb_transporters = len(self.transporters)
        self.nb_reversibles = len(self._reversibles_l)
        self.nb_nullreacs = len(self.nullreacs)
        

    def print_stats(self):
        """
        Prints the caracteristics of the model, such as the number of reactions and metabolites.
            
        """
        if hasattr(self, 'dual_mcs'): print("-- Dual network --")
        print(self.nb_metint, "internal metabolites")
        print(self.nb_metext, "external metabolites")
        print(self.nb_reactions, "total reactions")
        print(self.nb_transporters, "exchange reactions")
        print(self.nb_reversibles, "reversible reactions")
        
            
    def fill_smatrix(self):
        """
        Fills a stoichiometry matrix with the network data
        Stores all matrix related fields in the data structure
        
        Initializes:
            matrix: stoichiometry matrix
            ext_matrix: full stoichiometry matrix
            idx2reac, reac2idx: dicts to help reaction name/index association
            idx2imet, imet2idx: dicts to help metabolite name/index association
            nb_metint, nb_reactions: dimensions of the matrix
            nb_metabolites: alias for nb_metint, nb of internal metabolites
            ext_stoichiometry: sublist of stoichiometry only containing 
                               external metabolite stoichiometry tuples          
        """
        
        # Ensures self.metint is a list for performance
        self.metint = sorted(list(self.metint))

        # Matrix dimensions
        self.nb_reactions = len(self.reactions)
        self.nb_metint = len(self.metint)
        self.nb_metabolites = self.nb_metint
        
        # Additional dicts to help with direct association
        self.idx2met = {i: m for i, m in enumerate(self.metabolites)}
        self.idx2imet = {i: m for i, m in enumerate(self.metint)} 
        self.idx2reac = {i: r for i, r in enumerate(self.reactions)}  
        self.met2idx = {m: i for i, m in enumerate(self.metabolites)} 
        self.imet2idx = {m: i for i, m in enumerate(self.metint)}
        self.reac2idx = {r: i for i, r in enumerate(self.reactions)}
        
        # Stoichiometry matrix with nb_metint lines and nb_reactions columns
        self.matrix = np.zeros((self.nb_metint, self.nb_reactions), dtype=float)
        self.ext_matrix = np.zeros((len(self.metabolites), self.nb_reactions), dtype=float)
        self.ext_stoichiometry = []
        self.transporters = getattr(self, 'transporters', [])
        
        for tupl in self.stoichiometry:
            self.ext_matrix[self.met2idx[tupl[0]], self.reac2idx[tupl[1]]] = tupl[2]
            if tupl[0] in self.imet2idx:
                self.matrix[self.imet2idx[tupl[0]], self.reac2idx[tupl[1]]] = tupl[2]
            else:
                self.ext_stoichiometry.append(tupl)
                self.transporters.append(tupl[1])
        

    def __init__(self, metext, metabolites, reactions, stoichiometry, reversibles, nullreacs=[]):
        """
        Creates a MetaNetwork object from the parameters
        Stores all parameters in the data structure
        
        Params:
            metext: list of external metabolites names
            metabolites: list of metabolite names
            reactions: list of reactions names
            stoichiometry: list of stoichiometry tuples
                          (metabolite name, reaction name, value)
            reversibles: list of reversible reactions
            nullreacs: list of excluded exchange reactions
        
        Initializes:
            metint: set of internal metabolites names
            irreversibles: set of irreversible reactions names
            matrix: stoichiometry matrix
            transporters: list of transport reactions
         
        Calls method fill_smatrix and also initializes all matrix related fields
            
        """
        # Lists
        self.metabolites = sorted(list(metabolites))
        self.reactions = sorted(list(reactions))
        self.nullreacs = sorted(list(nullreacs))
        # Sets
        self.metext = set(metext)
        self.reversibles = set(reversibles)
        self.metint = set(self.metabolites).difference(self.metext)
        self.irreversibles = set(self.reactions).difference(self.reversibles)
        # Hidden utility structures
        self._metext_l = sorted(list(metext))
        self._metint_s = self.metint
        self._reversibles_l = sorted(list(reversibles))
        self._metabolites_s = set(metabolites)
        self._reactions_s = set(reactions)
        # Stoichiometry
        self.stoichiometry = stoichiometry
        self.fill_smatrix()
        # Transporters and extra fields
        self.transporters = sorted(list(set(self.transporters)))
        self.extras = getattr(self, 'extras', {})
        self.extras['transporter'] = self.transporters
        # Data structures check
        self.compute_stats()
        self.structures_check()
        self.null_reactions_check()
        self.only_one_stoch_check()
        # Print stats
        self.print_stats()

