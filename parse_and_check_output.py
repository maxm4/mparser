# -*- coding: utf-8 -*-

"""
parse_and_check_output.py
More advanced version of parse_output.py based on aspsmatrix's matrix
Uses kernel rank to check the output, verify if the flux modes are EFMs

"""

import numpy as np
import scipy.linalg
import numpy as np
import numpy.linalg
import pandas
import re, sys
from aspsmatrix import AspSmatrix, metatool_file, sbml_file, file_collection
from datetime import timedelta
from copy import deepcopy
from functools import reduce
from scipy.optimize import linprog
from argparse import ArgumentParser
from mparser import parsehelper as ph
import warnings
import pickle


SBML = 1
MTTL = -1
PICKLE = 0


def enumerate_list(lst):
    """
    Simple method enumerating a list and printing it to sysout

    Params: 
        lst: Python list
    """
    for i, elem in enumerate(lst):
        print(i+1, ":", elem)
        


def is_kernel_dimension_one(matrix):
    """
    Checks if the kernel of a matrix is of dimension one
    Uses SVD to calculate the rank of the matrix

    Params: 
        matrix: stoichiometric matrix

    Returns: 
        boolean: True if the kernel is of dimension one, else False
    """
    # utilise SVD à la place de Gaussian Elimination (Klamt, 2005)
    rank = numpy.linalg.matrix_rank(matrix)
    # print("forme matrice", matrix.shape, "rang", rank)
    if (rank == matrix.shape[1] - 1):  # nb columns aka nb reactions
        return True
    return False



def support_as_boolean(support, reacidx):
    """
    Converts support to a boolean table matching the network file matrix

    Params: 
        support: set of string names, support of the reaction, from clingoLP output
        reacidx: hash map associating reaction names to indices, from network file

    Returns: 
        booltable: table of booleans:
                   True if the reaction is in the support, else False
    """
    booltable = []
    itersup = support.copy()
    for val in reacidx: # for each reaction name
        match = val in itersup # boolean
        booltable.append(match)
        if match:
            itersup.remove(val)
    if itersup: # if support set is not empty
        raise ValueError("Reactions do not match between clingoLP and network file")
    return booltable



def get_efms_from_solutions(sols, matrix, reacidx):
    """
    Gets all EFMs from a list of solutions
    Uses the stoichiometric matrix
    Tests for each solution if the kernel of the submatrix is dimension one
    The submatrix is the stoichiometric matrix minus the inactive reactions columns 
    
    Params: 
        sols: Python List of solutions
        matrix: NumPy Array, stoichiometric matrix
        reacidx: Python Dict of names corresponding to matrix indices
        
    Returns:
        efms: Python List of EFMs
    """    
    efms = []
    for support in sols:
        if is_efm(support, matrix):
            efms.append(support)
    return efms 


def is_efm(support, matrix):
    """
    Checks if a support reaction names set is an elementary flux mode
    
    Params:
        fluxv: a flux vector returned by clingo[LP], as pandas dataframe
        matrix: original matrix, to check flux vectors minimality
                AspSmatrix format
    
    Returns:
        status: boolean, true if it is an efm else false
    """
    boolsupp = support_as_boolean(support, matrix.reactions_index())
    submatrix = matrix.matrix[:,boolsupp]
    return is_kernel_dimension_one(submatrix)


def str_is_efm(fluxv, matrix):
    """
    Checks if a flux vector is an elementary flux mode
    
    Params:
        fluxv: a flux vector returned by clingo[LP], as pandas dataframe
        matrix: original matrix, to check flux vectors minimality
                AspSmatrix format (optional)
    
    Returns:
        status: "unknown" if matrix is None, 
                "true" if fluxv is an elementary mode, else "false"
    """
    if fluxv.empty:
        raise ValueError("Empty flux vector dataframe")
    if matrix is None:
        return "unknown"
    else:
        support = list(map(ph.smart_remove_quotes, fluxv.index.values))
        support = set(list(map(ph.smart_remove_revs, support)))
        return str(is_efm(support, matrix)).lower()    

def get_flux_vectors(fname):
    """
    Reads the Clingo[LP] output and gets the flux vector for each solution
    TODO: Not use this function and get answers with clingo API instead
    
    Params: 
        fname: input file name, output from Clingo[LP]

    Returns:
        fluxes: list of flux vectors
        names: list of reaction or enzyme subset names
    """      
    f = open(fname)
    lines = f.readlines()
    fluxes = []
    names = []
    # uses the fact that the min function is 0 and always returns 0
    # after beginning brace index
    bb_idx = len('(0.0, {')
    # this filtering is a bit too specific to be reused...
    lines = filter(lambda w: w.startswith("(0.0,"), lines)
    spl = None
    length = -1
    for i, l in enumerate(lines):
        # end brace index
        eb_idx = l.rindex('})')
        # remove beginning and ending braces
        l = l[bb_idx:eb_idx]
        spl = l.split(",")
        try:
            ls = list(map(lambda w: float(w.split(":")[1]), spl))
            assert(length == -1 or length == len(spl))
            fluxes.append(np.array(ls))
        except (IndexError, AssertionError):
            warnings.warn(f'formatting defect on answer {i+1}')
        length = len(spl)

    if spl is None:
        raise ValueError() # Clingo LP output format error
    names = spl
    
    def get_word_in_quotes(e):
        res = re.compile(r'flux\(("?.*)"?\)').search(e)
        return "None" if res is None else res.group(1)
    
    names = list(map(get_word_in_quotes, names))

    f.close()
    return fluxes, names



def write_flux_vectors(inname, outname, checkfile, pklfile, mcfmfile, format, jsonfile, aspfile):
    """
    Writes the flux vectors retrieved from a clingo[LP] output
    
    Params: 
        inname: input file name, output from Clingo[LP]
        outname: file containing all elementary flux modes in text format
        checkfile: original network file, to check flux vectors minimality (optional)
        pklfile: additional file to store flux modes in pickle format (optional)
        mcfmfile: additional file to store non elementary flux modes indices in pickle format (optional)
        format: constant indicating if network file is in SBML format, METATOOL format or other
        jsonfile: additional file to store the last elementary flux mode in json format (optional)
        aspfile: additional file to store answer sets as an Answer Set Programming program
    """      
    fluxes, names = get_flux_vectors(inname)
    fvs = np.vstack(fluxes)
    pfvs = pandas.DataFrame(data=fvs, columns=names)
    pfvs = pfvs.reindex(sorted(pfvs.columns), axis=1)
    
    matrix = None
    if checkfile is not None:
        try:
            #format_file = sbml_file if sbml else metatool_file
            if format == SBML:
                matrix = AspSmatrix(sbml_file(checkfile))
            elif format == MTTL:
                matrix = AspSmatrix(metatool_file(checkfile))
            else:
                with open(checkfile, 'rb') as f:
                    data = pickle.load(f)
                    matrix = AspSmatrix(data, from_files=False)
        except IOError:
            print('File not found')
            
    if pklfile is not None:
        with open(pklfile, 'wb') as pkl:
            pickle.dump(pfvs, pkl)

    is_efm_ss = []
    fluxv = None
    with open(outname, "w") as f:
        for i in range(len(fvs)):
            fluxv = pfvs.loc[i, pfvs.loc[i]>0]
            f.write(str(fluxv))
            is_efm_str = str_is_efm(fluxv, matrix)
            is_efm_ss.append(is_efm_str)
            f.write(', elementary mode: ' + is_efm_str)
            f.write('\n\n')

    if mcfmfile is not None:
        is_efm_false = [i for i, x in enumerate(is_efm_ss) if x == 'false']
        with open(mcfmfile, 'wb') as pkl:
            pickle.dump(is_efm_false, pkl)

    if jsonfile is not None:
        with open(jsonfile, 'w') as j:
            fluxv.to_json(j)

    if aspfile is not None:
        raise NotImplementedError
  
if __name__== "__main__":
    parser = ArgumentParser()
    parser.add_argument('infile', metavar='clingoLP.file', help='Input file name')
    parser.add_argument('outfile', metavar='output.file', help='Output file name')
    parser.add_argument('--check', metavar='network.file', help='Check flux modes')
    parser.add_argument('--pickle', metavar='pickle.file', help='Store flux modes with pickle')
    parser.add_argument('--mcfm', metavar='mcfm.file', help='Store non elementary modes with pickle')
    parser.add_argument('--sbml', action='store_true', help='If network file is in SBML format')
    parser.add_argument('--metatool', action='store_true', help='If network file is in Metatool format')
    parser.add_argument('--json', metavar='json.file', help='Store flux modes in json')
    parser.add_argument('--solutions', metavar='asp.file', help='Store answer sets as an ASP program')
    opts = parser.parse_args()
    format = SBML if opts.sbml else MTTL if opts.metatool else PICKLE
    write_flux_vectors(opts.infile, opts.outfile, opts.check, opts.pickle, opts.mcfm, format, opts.json, opts.solutions)
