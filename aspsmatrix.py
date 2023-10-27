# -*- coding: utf-8 -*-

"""
aspsmatrix.py

Module for working with ASP and stoichiometric matrices
Reads stoichiometric matrices in ASP format
Computes enzyme subsets and other useful information

WARNING: Enzyme subsets don't work with decimal stoichiometry for now

"""

import numpy as np
import scipy.linalg
import sys, pickle 
import sympy, warnings
from pathlib import Path
from thermo_csv import read_csv_values
from mparser import InputFormatType, MetaNetwork
from mparser import read_network, format_asp, write_asp, append_asp
from mparser import parsehelper as ph
from argparse import ArgumentParser
from copy import deepcopy
from functools import reduce
import numpy.linalg
from sympy import Matrix as sympy_matrix
from sympy import sympify
from scipy.optimize import linprog
from fractions import Fraction


metatool_file = lambda f: file_collection(f, mtb_format=InputFormatType.TXT)
sbml_file = lambda f: file_collection(f, mtb_format=InputFormatType.SBML)


def get_names_collections(answers, quotes=True):
    """
    Gets collections of string names from a set of names
    
    Params: 
        answers: set of names
        quotes: boolean, indicating whether names are in quotes
        
    Returns: 
        collections: dict containing various collections
    """        
    collections = {}
    collections['map'] = {}
    collections['map_nq'] = {}
    collections['map_inv'] = {}
    collections['map_inv_nq'] = {}
    collections['set'] = set()
    collections['set_nq'] = set()
    
    remove_quotes = lambda w : (w, w[1:-1]) # name_nq string with no quotes
    add_quotes = lambda w : ('"{}"'.format(w), w) # name string with added quotes
    handle_quotes = remove_quotes if quotes else add_quotes
    
    for i, value in enumerate(sorted(answers)):
        name = value # gets the name string, presumably the only element from the tuple
        name, name_nq = handle_quotes(name)
        collections['map'][name] = i
        collections['map_nq'][name_nq] = i
        collections['map_inv'][i] = name
        collections['map_inv_nq'][i] = name_nq
        collections['set'].add(name)
        collections['set_nq'].add(name_nq)
        
    return collections


def file_collection(network_file, mtb_format=InputFormatType.SBML,
                    constraints_file=None, keq_thermo_file=None,
                    extc_thermo_file=None):
    """
    Creates a collection of file names to be read and written
    
    Params: 
        network_file: input metabolism network file name
        mtb_format: metabolism network format (default: SBML)
        constraints_file: file containing arbitrary ASP constraints
        keq_thermo_file: file containing reaction equilibrium constants
        extc_thermo_file: file containing external metabolite concentrations
        
    Returns: 
        collection: file collection dict containing all these parameters as well as
                    auto-named ASP output files for network, constraints and subsets
    """            
    collection = {}
    
    collection['network_file'] = network_file
    collection['mtb_format'] = mtb_format
    collection['constraints_file'] = constraints_file
    collection['keq_thermo_file'] = keq_thermo_file
    collection['extc_thermo_file'] = extc_thermo_file
    
    collection['name'] = Path(network_file).stem
    collection['asp_result_file'] = collection['name'] + ".compressed.lp4"
    collection['asp_constraints_file'] = collection['name'] + ".constraints.lp4"
    collection['asp_subsets_file'] = collection['name'] + ".subsets.lp4"
    collection['asp_left_nullspace'] = collection['name'] + ".left_nullspace.lp4"
    collection['pkl_left_nullspace'] = collection['name'] + ".left_nullspace.pkl"
    collection['pkl_aspsmatrix'] = collection['name'] + ".aspsmatrix.pkl"

    return collection


class AspSmatrix:
    """
    
    AspSmatrix class
    Class for handling stoichiometric matrices in ASP
    
    """ 
    
    def __init__(self, files, no_external=False, from_files=True):
        """
        Creates an AspSmatrix object from reading an ASP file
        
        Params:
            files: file collection dict containing the following parameters:
                network_file: input metabolism network file name
                mtb_format: metabolism network format (default: txt)
                constraints_file: file containing arbitrary ASP constraints
                keq_thermo_file: file containing reaction equilibrium constants
                extc_thermo_file: file containing external metabolite concentrations
            no_external: boolean indicating the presence of external stoichiometry
                         if set to true external stoichiometry is considered internal
        """ 
        self.from_files = from_files

        if self.from_files:
            self._network_file = files['network_file']
            self._mtb_format = files['mtb_format']
            self._constraints_file = files['constraints_file']
            self._keq_thermo_file = files['keq_thermo_file']
            self._extc_thermo_file = files['extc_thermo_file']
        else:
            #assert(type(files) == MetaNetwork)
            self._mtb_network = files             

        self._quotes = False
        self.nullspace = False
        self._read_network_(no_external)
        self._read_constraints_()
        self.reduced_network = None

        if self.nullspace:
            print("Computing stoichiometric kernel")
            print(self._ns_symb)
            

    def _read_network_(self, no_external):
        """
        Reads a metabolic network and converts it to a stoichiometric matrix
        
        Params:
            no_external: boolean indicating the presence of external stoichiometry
                         if set to true external stoichiometry is considered internal
        
        Uses: 
            network_file: input file name computed to get the stoichiometric matrix
            mtb_format: metabolism network format (default: txt)
            
        Initializes:
            all_metabolites: set of all metabolites, including external metabolites
            reactions: dict containing collections of reactions names and indices
            metabolites: dict containing collections of metabolites names and indices
                         external metabolites are excluded from this dict and the matrix
            reversibles: set of reversible reactions
            is_reversible: list of boolean values
            nb_reactions: number of reactions, number of columns of the matrix
            nb_metabolites: number of metabolites, number of lines of the matrix
            stoichiometry: full stoichiometry graph, including external metabolites
                           list of (metabolite name, reaction name, value) tuples
            ext_stoichiometry: external metabolites stoichiometry graph
            matrix: stoichiometry matrix, NumPy Array of dtype int
            ns_symb: kernel matrix, NumPy Array of SymPy symbols
            ns_float: kernel matrix, NumPy Array of decimals
        """ 
        if self.from_files:
            print("Reading file", self._network_file)
            netw = read_network(self._network_file, self._mtb_format)
        else: 
            netw = self._mtb_network
        
        self._all_metabolites = netw.metabolites # external metabolites are included 
        self._metext = netw.metext # external metabolites   
        
        netw.metint = netw.metabolites if no_external else netw.metint        
        self._metabolites = get_names_collections(netw.metint, self._quotes)  
        self._reactions = get_names_collections(netw.reactions, self._quotes)
        
        self.nb_reactions = len(netw.reactions)
        self.nb_metabolites = len(netw.metint)
        
        self._reversibles = netw.reversibles
        self._stoichiometry = netw.stoichiometry
            
        metaidx = self.metabolites_index()   
        reacidx = self.reactions_index()   
        invreacidx = self.reactions_index(reverse=True)  
        
        self._is_reversible = [invreacidx[x] in self._reversibles 
                               for x in range(self.nb_reactions)]
            
        # stoichiometry matrix of type int with nb_metabolites lines and nb_reactions columns
        self.matrix = np.zeros((self.nb_metabolites, self.nb_reactions), dtype=float)
        self._ext_stoichiometry = []
        
        for tupl in self._stoichiometry: # external metabolites may have stoichiometry
            if tupl[0] in metaidx:  # checks if this is an internal metabolite
                self.matrix[metaidx[tupl[0]], reacidx[tupl[1]]] = tupl[2]
            else:
                self._ext_stoichiometry.append(tupl) # used for thermo constraints
          
        if self.nb_reactions > 50:
            warnings.warn("large model; skipping nullspace computation")
            return
           
        try:
            self.nullspace_computation()
        except ValueError as error: 
            warnings.warn(str(error)) # ValueError
        except OverflowError as error: 
            warnings.warn(str(error)) # OverflowError


    def nullspace_computation(self, skip_ns_int=False):
        """
        Performs the integer nullspace computation if skipped
        Computationally expensive

        Throws: ValueError

        Params:
            skip_ns_int: If True, skips integer nullspace computation

        Initializes:
            ns_symb, ns_float, ns_int: Nullspace in various types
            nullspace: Sets the flag to true
        """
        if self.nullspace:
            return       
        self._ns_int = None
        self._ns_symb, self._ns_float = self.compute_nullspace(self.matrix) 
        if not skip_ns_int:
            self._ns_int = self.compute_integer_nullspace(self._ns_symb)
        self.nullspace = True
                   
    
    def _read_constraints_(self):
        """
        Reads additional constraints for the metabolic network
        
        Uses: 
            constraints_file: file containing arbitrary ASP constraints
            keq_thermo_file: file containing reaction equilibrium constants
            extc_thermo_file: file containing external metabolite concentrations
            
        Initializes:
            names_keq: list of reaction names
            values_keq: list of apparent equilibrium constants
                        adapted for positive flux values            
        Note:
            Both initialized lists have the same indexing
        """ 
        if not self.from_files:
            return

        if self._constraints_file is not None:
            print("With constraints", self._constraints_file)
            
        if self._keq_thermo_file is not None:
            print("With equilibrium constants", self._keq_thermo_file)
            reac_names, keq = read_csv_values(self._keq_thermo_file)
            
            if self._extc_thermo_file is not None:
                print("And external concentrations", self._extc_thermo_file)
                metext_names, extc = read_csv_values(self._extc_thermo_file)
                prod = np.ones(len(keq)) # product of metext concentrations
                for tupl in self._ext_stoichiometry:
                    meta, reaction, stoch = tupl # gets stoichiometry
                    m = ph.indexl(metext_names, meta)
                    if m != -1: # if metabolite in meta_names
                        r = ph.indexl(reac_names, reaction)
                        if r != -1: # if reaction in meta_names
                            if extc[m] != 0: 
                                prod[r] *= (extc[m] ** stoch)
                        else:
                           print("Unmatched reaction name", reaction) 
                    else:
                        print("Unmatched metabolite name", meta)
                        
                keq /= prod
                
            else:
                print("No external concentrations file found")
                print("Equilibrium constants will be used directly")
                
            # neparian logarithm
            self._values_keq = np.log(keq)
            self._names_keq = reac_names    
    
    
    def _index_str_(self, quotes, reverse, as_set):
        """
        Gets the name of the index collection depending on parameters
        
        Params: 
            quotes: boolean determining whether dict entries have quotes or not
                    if not given, determines the appropriate choice with self._quotes
            reverse: boolean determining if the returned map is mapping indices to
                     reactions names (if True) or reactions names to indices (if False)
            as_set: boolean determining if the index should be a dict (False) or a set (True)
            
        Returns: 
            name: name of the corresponding collection
        """ 
        if quotes == None:
            quotes = self._quotes
        if quotes == True:
            if as_set == True:
                return 'set'
            if reverse == True:
                return 'map_inv' 
            else:
                return 'map'
        if as_set == True:
            return 'set_nq'
        if reverse == True:
            return 'map_inv_nq' 
        else:
            return 'map_nq'       
    
    
    def metabolites_index(self, quotes=None, reverse=False, as_set=False):
        """
        Gets a dict associating metabolites names to their index in the stoichiometric matrix
        
        Params: 
            quotes: boolean determining whether dict entries have quotes or not
                    if not given, determines the appropriate choice with self._quotes
            reverse: boolean determining if the returned map is mapping indices to
                     reactions names (if True) or reactions names to indices (if False)
            as_set: boolean determining if the index should be a dict (False) or a set (True)
            
        Returns: 
            metaidx: dict associating metabolites names to indices
                     set of metabolites names if as_set is True
        """  
        return self._metabolites[self._index_str_(quotes, reverse, as_set)]
    
    
    def reactions_index(self, quotes=None, reverse=False, as_set=False):
        """
        Gets a dict associating reactions names to their index in the stoichiometric matrix
        
        Params: 
            quotes: boolean determining whether dict entries have quotes or not
                    if not given, determines the appropriate choice with self._quotes
            reverse: boolean determining if the returned map is mapping indices to
                     reactions names (if True) or reactions names to indices (if False)
            as_set: boolean determining if the index should be a dict (False) or a set (True)
            
        Returns: 
            reacidx: dict associating reactions names to indices
                     set of reactions names if as_set is True
        """
        return self._reactions[self._index_str_(quotes, reverse, as_set)] 


    def compute_integer_nullspace(self, symbmat):
        """
        Converts SymPy rationals and fractions to integers
        Also simplifies each integer row by dividing by their gcd
        In a separate function because the operation costs computation power
        
        Params: 
            symbmat: kernel matrix, NumPy Array of SymPy rational symbols 
        
        Returns:
            ns_int: kernel matrix, NumPy Array of integers
        """
        matrix = np.zeros(symbmat.shape, dtype=np.int64)
        lcm = 1
        for i, row in enumerate(symbmat):
            denom_row = list(map(lambda e: e.as_numer_denom()[1], row))
            lcm = np.lcm(lcm, np.lcm.reduce(denom_row)) # least common multiple
        for i, row in enumerate(symbmat):
            matrix[i] = row * lcm # convert fractions row to an int row
        return matrix
    
    
    def compute_old_integer_nullspace(self, symbmat):
        """
        Converts SymPy rationals and fractions to integers
        Also simplifies each integer row by dividing by their gcd
        In a separate function because the operation costs computation power
        
        Params: 
            symbmat: kernel matrix, NumPy Array of SymPy rational symbols 
        
        Returns:
            ns_int: kernel matrix, NumPy Array of integers
        """
        matrix = np.zeros(symbmat.shape, dtype=int)
        for i, row in enumerate(symbmat):
            denom_row = list(map(lambda e: e.as_numer_denom()[1], row))
            lcm = np.lcm.reduce(denom_row) # least common multiple
            matrix[i] = row * lcm # convert fractions row to an int row
        return matrix
    
    
    def compute_nullspace(self, matrix):
        """
        Computes a kernel matrix from the nullspace of a given matrix
        
        Params: 
            matrix: a given matrix, NumPy Array, must be of dtype int
        
        Returns:
            ns_symb: kernel matrix, NumPy Array of SymPy symbols
            ns_float: kernel matrix, NumPy Array of decimals
        """
        f = lambda d: Fraction.from_float(d).limit_denominator(10000)
        f = np.vectorize(f)
        matrix = f(matrix)
        smatrix = sympy.Matrix(matrix)
        nullspace = smatrix.nullspace(simplify=False)
        if not nullspace:
            raise ValueError("Empty nullspace")
        ns_symb = np.hstack(sympy.list2numpy(nullspace))
        ns_float = np.array(ns_symb, dtype=float)
        return ns_symb, ns_float
    
    
    def enzyme_subsets_vector(self):
        """
        Computes enzyme subsets from the matrix and gives the results as a vector
            
        Returns: 
            array: column vector indicating for each row to which enzyme subset they belong
            nb_subsets: final number of subsets
        """
        if not self.nullspace:
            raise ValueError("Requires stoichiometric kernel")        
        ker = self._ns_int
        rset = {} # dict mapping rows to subsets number
        array = [0] * len(ker)
        current_subset = 1
        for i, row in enumerate(ker):
            gcd = np.gcd.reduce(row)
            if gcd != 0:
                row = row / gcd # normalize column to compare it with others
                htr = str(row) # uses string representation of row to compare
                if htr in rset:
                    array[i] = rset[htr]
                else:
                    rset[htr] = current_subset
                    array[i] = current_subset
                    current_subset += 1
        return array, current_subset       
    
    
    def get_alpha_factor(self, row1, row2):
        """
        Gets the factor alpha by which the two kernel lines are colinear
        If the two lines are not colinear, returns None
            
        Returns: 
            alpha: SymPy rational, colinearity factor, else None
        """
        nrow1 = row1[np.nonzero(row1)]
        nrow2 = row2[np.nonzero(row2)]
        div = nrow1 / nrow2
        elem = div[0]
        if not np.all(div == elem):
            raise ValueError("The two given rows are not colinear")
        return elem
    
    
    def verify_enzyme_subsets(self, subsets_vector, nb_subsets):
        """
        Verifies if enzyme subsets are valid and calculates factor coefficients
        
        Params:
            subsets_vector: list indicating for each row the enzyme subset
            nb_subsets: number of subsets, useful for declaring the list of subsets
        
        Returns: 
            to_delete: set of enzyme subsets to be deleted
            es_reversible: list of booleans indicating that subsets are reversibles if True
            coeff_subsets: enzyme subsets as a list of lists [factor, index]
                           where factor is a multiplier for the stoichiometry
            
        """
        if not self.nullspace:
            raise ValueError("Requires stoichiometric kernel")

        # search alpha factor coefficients            
        firstin = [-1 for i in range(nb_subsets)] # for each subset
        es_reversible = [True for i in range(nb_subsets)] # for each subset
        coeff_subsets = [[] for i in range(nb_subsets)]
        to_delete = {0} # enzyme subset 0 corresponds to zero rows and are always deleted
        for r, subset in enumerate(subsets_vector):
            if subset in to_delete:
                continue
            if firstin[subset] == -1:
                firstin[subset] = r
                coeff_subsets[subset].append([1, r])
                es_reversible[subset] = self._is_reversible[r]
            else:
                first = firstin[subset]
                alpha = self.get_alpha_factor(self._ns_symb[r], self._ns_symb[first])
                coeff_subsets[subset].append([alpha, r])
                es_reversible[subset] = es_reversible[subset] and self._is_reversible[r]
                # no flux if both reactions are irreversible with negative alpha
                if alpha < 0 and not (self._is_reversible[r] or self._is_reversible[first]): 
                    to_delete.add(subset)
        
        # may not be great time complexity wise
        for subset in coeff_subsets:
            if len(subset) > 1:
                denom_row = list(map(lambda e: sympify(e[0]).as_numer_denom()[1], subset))
                lcm = np.lcm.reduce(denom_row) # least common multiple
                for tupl in subset:
                    tupl[0] = int(tupl[0] * lcm)
        
        return to_delete, es_reversible, coeff_subsets
    
    
    def make_reduced_network(self):
        """
        Computes a reduced network given the enzyme subsets
        Combines stoichiometry of reactions in the same subset
        
        Uses:
            subsets_vector: list indicating for each row the enzyme subset
            to_delete: set of enzyme subsets to be deleted
            es_reversible: list of booleans indicating that subsets are reversibles if True
            coeff_subsets: enzyme subsets as a list of lists [factor, index]
        
        Returns: 
            network: reduced network, MetaNetwork object
            
        Initializes:
            reduced_reactions: set of reduced reactions names
            reduced_reversibles: set of reduced reversible reactions names
            subsets_vector: list indicating for each row the enzyme subset
            all_metabolites: set of all metabolites, including external metabolites
            nb_subsets: number of subsets
            coeff_subsets: enzyme subsets as a list of lists [factor, index]
            metabolites: dict containing collections of metabolites names and indices
                         external metabolites are excluded from this dict and the matrix
            reduced_matrix: reduced stoichiometry matrix, NumPy Array of dtype int
            reduced_stoichiometry: reduced stoichiometry graph, without external metabolites
                                   list of (metabolite name, reaction name, value) tuples
        """
        if not self.nullspace:
            raise ValueError("Requires stoichiometric kernel")
        
        self.subsets_vector, total_nb_subsets = self.enzyme_subsets_vector() 
        verified_subsets = self.verify_enzyme_subsets(self.subsets_vector, total_nb_subsets)
        to_delete, es_reversible, self._coeff_subsets = verified_subsets
        
        self._reduced_reactions = set()
        self._reduced_reversibles = set()
        for sub in range(total_nb_subsets):
            if sub not in to_delete:
                self._reduced_reactions.add(str(sub))
                if es_reversible[sub]:
                    self._reduced_reversibles.add(str(sub))
        self.nb_subsets = len(self._reduced_reactions)
                    
        invmetaidx = self.metabolites_index(reverse=True)
        invsubsetidx = {}
        
        j = 0
        self.reduced_matrix = np.zeros((self.nb_metabolites, self.nb_subsets), dtype=float)
        for i, sub in enumerate(self._coeff_subsets):
            if i in to_delete:
                continue
            for elem in sub:
                coeff, idx = elem
                
                self.reduced_matrix[:, j] += coeff * self.matrix[:, idx]
            invsubsetidx[j] = str(i)
            j += 1 # i indexes all subsets while j indexes only valid subsets
                        
        self._reduced_stoichiometry = []
        nonzero = self.reduced_matrix.nonzero() # is it really better complexity-wise?
        for i in range(len(nonzero[0])) : # external metabolites may have stoichiometry
            x = nonzero[0][i]
            y = nonzero[1][i]
            meta = invmetaidx[x]
            reaction = invsubsetidx[y]
            self._reduced_stoichiometry.append((meta, reaction, self.reduced_matrix[x][y]))    
        
        metext = self._metext
        metabolites = self._all_metabolites
        reactions = self._reduced_reactions
        stoichiometry = self._reduced_stoichiometry
        reversibles = self._reduced_reversibles
        self.reduced_network = MetaNetwork(metext, metabolites, reactions,
                                           stoichiometry, reversibles)
        
            
    def feasibility_analysis(self):
        """
        Feasibility analysis using Linear Programming
        As described in the following PhD Thesis by Marco Terzer:
        "Large Scale Methods to Enumerate ExtremeRays and Elementary Modes"
        Warning: NOT WORKING YET
            
        Returns: 
            zero_flux: set of reactions with zero flux value, may be removed
            îrreversibles_reversibles: set of reversibles reactions that may be irreversible
            essential_reactions: set of reactions that are essential
        """
        zero_flux = set()
        irreversibles_reversibles = set()
        essential_reactions = set()
        opt = numpy.ones((self.nb_reactions,), dtype=int)
        cstr = numpy.zeros((self.nb_metabolites,), dtype=int)
        rbounds = [() for i in range(self.nb_reactions)]
        for i in range(self.nb_reactions):
            if self._is_reversible[i]:
                rbounds[i] = (None, None)
            else:
                rbounds[i] = (0, None)
        minres = linprog(opt,  A_ub=self.matrix, b_ub=cstr, bounds=rbounds)
        maxres = linprog(-opt, A_ub=self.matrix, b_ub=cstr, bounds=rbounds)
        if minres.success != True or maxres.success != True:
            raise ValueError("Linear programming failed to find an optimal solution")
        for i in range(self.nb_reactions):
            if minres.x[i] == 0 and maxres.x[i] == 0:
                zero_flux.add(i)
            if self._is_reversible[i] and (minres.x[i] == 0 or maxres.x[i] == 0):
                irreversibles_reversibles.add(i)
            if minres.x[i] != 0 and maxres.x[i] != 0 and np.sign(minres.x[i]) == np.sign(maxres.x[i]):
                essential_reactions.add(i)
        
        return zero_flux, irreversibles_reversibles, essential_reactions
    

    def _constraints_to_asp_(self):
        """
        Creates ASP rules for query constraints [TODO]
            
        Returns: 
            constraints: string containing the constraints
        """
        if self._constraints_file is None:
            return ""
        constraints = ""
        return constraints
    

    def _thermo_to_asp_(self):
        """
        Creates ASP rules for thermodynamic constraints
        
        Uses:
            names: list of reaction names 
            values: list of log apparent equilibrium constants
            reactions: set of reaction names
            reversibles: set of reversible reaction names
            
        Returns: 
            constraints: string containing the constraints
        """
        if self._keq_thermo_file is None: # then names and values aren't defined
            return ""
        constraints = ""
        names = self._names_keq
        values = self._values_keq
        withquotes = lambda w: f'"{w}"'
        for i in range(len(names)):
            name = names[i]
            if name in self.reactions_index(as_set=True):
                svalue = str(values[i])
                qname = withquotes(name)
                value = withquotes(svalue)
                constraints += "keq(" + qname + "," + value + ").\n"
                if name in self._reversibles:
                    qname_rev = withquotes(name + "_rev")
                    value = withquotes(svalue)
                    constraints += "keq(" + qname_rev + "," + value + ").\n"
            else:
                print("Ignored reaction name", name)
        return constraints
    
    
    def write_asp_network(self, outfname):
        """
        Writes the network as ASP rules in the specified file
        
        Params:
            outfname: output file name
        """
        if self.reduced_network is None:
            self.make_reduced_network()
        format_asp(self.reduced_network, outfname)
        print("Generated file", outfname)    
            
        
    def write_asp_constraints(self, outfname):
        """
        Writes the constraints as ASP rules in the specified file
        If there is no additional constraints provided, writes a blank file
        
        Params:
            outfname: output file name
        """
        constraints = self._constraints_to_asp_() + self._thermo_to_asp_()
        if constraints:
            write_asp(constraints, outfname)
            print("Generated file", outfname)
        
        
    def convert_reduced_network_solutions(self, fluxes, names):
        """
        Converts the reduced network solutions to full network solutions
        
        Params:
            fluxes: list of flux vectors, solutions of the reduced network
            names: list of reaction or enzyme subset names
        
        Returns:
            solutions: list of flux vectors dicts, solutions of the full network
                       a flux is represented as a dict with key/value = names/fluxes
        
        Note:
            Returns None if called before make_reduced_network
        """
        if self.reduced_network is None:
            return None
        invreacidx = self.reactions_index(reverse=True)
        reversibles = list(map(lambda a: True if "_rev" in a else False, names))
        names = list(map(lambda a: int(a.split("_")[0]), names))
        finalfluxes = [{} for i in range(len(fluxes))] # table of dicts
        for i, flux in enumerate(fluxes):
            for j, val in enumerate(flux):
                if val != 0.0:
                    for duo in self._coeff_subsets[names[j]]:
                        coeff, reac = duo
                        if reversibles[j] :
                            coeff *= -1
                        name = invreacidx[reac]
                        finalfluxes[i][name] = coeff * val       
        return finalfluxes
        
    
def compute_neighbours(aspsm):
    """
    Compute a list of neighbours for each reaction of the aspsmatrix
    Considers reversible reactions too
    
    Params:
        aspsm: AspSmatrix object
        
    Returns: 
        neighbours: dict of sets of neighbours for each reaction
    """
    neighbours = {}
    ridex = aspsm.reactions_index()
    idexr = aspsm.reactions_index(reverse=True) 
    for reac in ridex.keys():
        neighbours[reac] = set()
        nz = list(aspsm.matrix[:,ridex[reac]].nonzero()[0])
        for meta in nz:
            for neighbour in aspsm.matrix[meta].nonzero()[0]:
                neighbours[reac].add(idexr[neighbour])
        if reac in neighbours[reac]: 
            neighbours[reac].remove(reac)
    return neighbours


def asp_enzyme_subset(nid, reaction):
    """
    Creates an enzyme subset ASP constraint to add to the ASP solver
    
    Params:
        nid: number representing the subset 
        reaction: reaction name
        
    Returns: 
        enzyme_subset: string containing the ASP constraint
    """
    return "enzyme_subset(" + str(nid) + "," +  reaction + ").\n"



def enzyme_subsets(aspsm, files, query=None):
    """
    Computes enzyme subsets
    
    Params: 
        aspsm: AspSmatrix structure
        files: file collection as defined in the method
        query:  additional file containing a query, else None
    """        
    aspsm.nullspace_computation()
    invreacidx = aspsm.reactions_index(reverse=True)               
    aspsm.make_reduced_network()
    # writes network into ASP
    aspsm.write_asp_network(files['asp_result_file'])
    aspsm.write_asp_constraints(files['asp_constraints_file'])
    subsets_vector = aspsm.subsets_vector
    #print("reacidx", reacidx)
    #print(subsets_vector)
    #print(aspsm.matrix)
    
    enzyme_subsets_str = "" 
    for i in range(0, len(subsets_vector)) :
        enzyme_subsets_str += asp_enzyme_subset(subsets_vector[i], invreacidx[i])
       
    out = open(files['asp_subsets_file'], "w")
    out.write(enzyme_subsets_str)
    out.close()
    
    print("Generated file", out.name)
    
    
    
def left_nullspace(aspsm, files):
    """
    Computes left nullspace
    
    Params: 
        aspsm: AspSmatrix structure
        files: file collection as defined in the method
    """   

    metaidx = aspsm.metabolites_index() 
    invmetaidx = aspsm.metabolites_index(reverse=True)    
    #print(invmetaidx) 
    
    lns_symb, lns_float = aspsm.compute_nullspace(aspsm.matrix.T)
    #print(lns_symb)
    lns_int = aspsm.compute_integer_nullspace(lns_symb)
    #print(lns_int)
    
    with open(files['pkl_left_nullspace'], 'wb') as out: 
        pickle.dump(lns_int, out)
        print("Generated file", out.name)
        
    asp = ""
    for idx, meta in invmetaidx.items(): # row
        for col, coeff in enumerate(lns_int[idx]): # column
            asp += f"leftkernel({ph.smart_add_quotes(meta)}, {col}, {coeff}).\n"
            
    asp += f"column(0..{len(lns_int[idx])-1}).\n"
        
    with open(files['asp_left_nullspace'], 'w') as out: 
        out.write(asp)
        print("Generated file", out.name)
   
"""     
    unconserved = {"AP", "BP", "CP"}    
    Y = set()
    for idx, meta in invmetaidx.items():
        if meta not in unconserved:
            continue
        if not np.any(lns_int[idx]):
            y = np.zeros(np.size(lns_int, 0))
            y[idx] = 1
            Y.add(y)
            print("met", meta)
        else:
            # MIP
            print(lns_int[idx])
            print("met", meta)
"""            
    
    
def asp_smatrix(opts):
    """
    Main program
    Reads a file, convert it to an AspSmatrix structure
    Saves it to a pickle file
    
    Params: 
        opts: parameters from ArgumentParsers, including:
            infile: input file name, for defining a file collection
            subsets: boolean indicating computation of enzyme subsets
            left_nullspace: boolean indicating computation of left nullspace
            
    """  
    format_collection = sbml_file if opts.sbml else metatool_file
    collec = format_collection(opts.infile)
    aspsm = AspSmatrix(collec, opts.no_external)
    #print(aspsm.matrix)
    
    if opts.subsets:
        enzyme_subsets(aspsm, collec, None)
        
    if opts.left_nullspace:
        left_nullspace(aspsm, collec)
        
    with open(collec['pkl_aspsmatrix'], 'wb') as out: 
        pickle.dump(aspsm, out)       
        print("Generated file", out.name)

    
  
if __name__== "__main__":
    parser = ArgumentParser()
    parser.add_argument('infile', metavar='network.file', help='Input file name')
    parser.add_argument('--subsets', action='store_true', help='Computes enzyme subsets')
    parser.add_argument('--left-nullspace', action='store_true', help='Computes left nullspace')
    parser.add_argument('--sbml', action='store_true', help='If network file is in SBML format')
    parser.add_argument('--no-external', action='store_true', help='If external stoichiometry is not wanted')
    opts = parser.parse_args()
    asp_smatrix(opts)


