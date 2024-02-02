#!/usr/bin/env python
# coding: utf-8

"""
compress_model.py

Compresses an input metabolic model in MetaNetwork format
Coded to be integrated into mparser cli tool

Authors: Maxime Mahout for all non-efmtool-link code
         Philip Schneider / CNAPy team for all efmtool-link code

"""

from argparse import ArgumentParser

import os, numpy, optlang
from typing import Dict
import tempfile
import cobra, cobra.util.array

import pprint 

import efmtool_link.efmtool_intern as efmtool_intern
import efmtool_link.efmtool_extern as efmtool_extern
import efmtool_link.efmtool4cobra as efmtool4cobra

from copy import deepcopy
import numpy as np
from .meta_network import MetaNetwork
from .network_to_sbml import format_sbml
from .metatool_to_network import txt_read_network

from pathlib import Path
import pandas, io, pickle

#############################################

def limit(word, max_letters):
    return word[:max_letters] + '...' if len(word) > max_letters else word

def numeric_naming_scheme(i, rxns, starting_index=1):
    return 'rsub_{}'.format(i+starting_index)

def starry_combination(i, rxns, starting_index=1):
    return limit('*'.join(rxns), 80)

naming_scheme = starry_combination # numeric_naming_scheme

#############################################

def flux_variability_analysis(model: cobra.Model, loopless=False, fraction_of_optimum=0.0,
                              processes=None, results_cache_dir: Path=None, fva_hash=None, print_func=print):
    """
    Performs FVA on a model, enhanced version of cobra's fva using hashmap to only run the operation if not already ran
    
    CODE FROM EFMTOOL LINK NOTEBOOK by CNAPy team
    """
    # all bounds in the model must be finite because the COBRApy FVA treats unbounded results as errors
    if results_cache_dir is not None:
        fva_hash.update(pickle.dumps((loopless, fraction_of_optimum, model.tolerance))) # integrate solver tolerances?
        if fraction_of_optimum > 0:
            fva_hash.update(pickle.dumps(model.reactions.list_attr("objective_coefficient")))
            fva_hash.update(model.objective_direction.encode())
        file_path = results_cache_dir / (model.id+"_FVA_"+fva_hash.hexdigest())
        fva_result = None
        if Path.exists(file_path):
            try:
                fva_result = pandas.read_pickle(file_path)
                print_func("Loaded FVA result from", str(file_path))
            except:
                print_func("Loading FVA result from", str(file_path), "failed, running FVA.")
        else:
            print_func("No cached result available, running FVA...")
        if fva_result is None:
            fva_result = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=None, loopless=loopless,
                                                             fraction_of_optimum=fraction_of_optimum,
                                                             pfba_factor=None, processes=processes)
            try:
                fva_result.to_pickle(file_path)
                print_func("Saved FVA result to ", str(file_path))
            except:
                print_func("Failed to write FVA result to ", str(file_path))
        return fva_result
    else:
        return cobra.flux_analysis.flux_variability_analysis(model, reaction_list=None, loopless=loopless,
                                                             fraction_of_optimum=fraction_of_optimum,
                                                             pfba_factor=None, processes=processes)


def fva_model(model, compressed_model):
    """
    Performs FVA on the model to be compressed to remove blocked reactions
    Constrains flux bounds of the model to their FVA values
    
    CODE FROM EFMTOOL LINK NOTEBOOK by CNAPy team
    """
    fva_tolerance=1e-9
    with model as fva: # can be skipped when a compressed model is available
        # when include_model_bounds=False modify bounds so that only reversibilites are used?
        # fva.solver = 'glpk_exact' # too slow for large models
        fva.tolerance = fva_tolerance
        fva.objective = model.problem.Objective(0.0)
        if fva.problem.__name__ == 'optlang.glpk_interface':
            # should emulate setting an optimality tolerance (which GLPK simplex does not have)
            fva.solver.configuration._smcp.meth = optlang.glpk_interface.GLP_DUAL
            fva.solver.configuration._smcp.tol_dj = fva_tolerance
        elif fva.problem.__name__ == 'optlang.coinor_cbc_interface':
            fva.solver.problem.opt_tol = fva_tolerance
        fva_res = flux_variability_analysis(fva, fraction_of_optimum=0.0, processes=1, 
                                                 results_cache_dir=None, fva_hash=None)
    for i in range(fva_res.values.shape[0]): # assumes the FVA results are oreduced_matrixered same as the model reactions
        if abs(fva_res.values[i, 0]) > fva_tolerance: # resolve with glpk_exact?
            compressed_model.reactions[i].lower_bound = fva_res.values[i, 0]
        else:
            compressed_model.reactions[i].lower_bound = 0
        if abs(fva_res.values[i, 1]) > fva_tolerance: # resolve with glpk_exact?
            compressed_model.reactions[i].upper_bound = fva_res.values[i, 1]
        else:
            compressed_model.reactions[i].upper_bound = 0
    return compressed_model
    
    
def get_reactions_equations(df, model):
    """
    Compute reaction equations for every reaction
    Uses StringIO for string writing and concatenation
    """
    with io.StringIO() as output:
        for column in df.columns:
            dfc = df[column]
            ftime=True
            print(column, ': ', end='', file=output, flush=True)
            for k, v in dfc[dfc < 0].items():
                if ftime:
                    ftime=False
                else:
                    print('+', end=' ', file=output, flush=True)
                print(-v, k, end=' ', file=output, flush=True)
            print('=' if model.reactions.get_by_id(column).reversibility else '=>', end=' ', file=output, flush=True)
            ftime=True
            for k, v in dfc[dfc > 0].items():
                if ftime:
                    ftime=False
                else:
                    print('+', end=' ', file=output, flush=True)
                print(v, k, end=' ', file=output, flush=True)
            print(file=output, flush=True)
        lines_r = output.getvalue()
    return(lines_r)
    
 
def make_compression_dict(model, original_model, subT, rxn_in_sub):
    """
    Makes correspondance dict between reaction and reaction subsets
    Also suggests bounds for ASP from compressed model FVA
    Bounds for reversibles (-LB, UB) are split: (0, UB), rev: (0, LB) 
    """
    ri = lambda i: original_model.reactions[i].id
    subname = lambda i, rxns: naming_scheme(i, rxns)
    subnamerev = lambda i, rxns: rev(subname(i, rxns))
    aq = lambda x: '"{}"'.format(x)
    aqrev = lambda x: '"{}_rev"'.format(x)
    rev = lambda x: '{}_rev'.format(x)
    dicto = {}
    dicto_keys = [] 

    for r, (ls, rc) in enumerate(zip(rxn_in_sub, model.reactions)):
        if len(ls) > 0:
            #print(ls, [(ri(rk), subT[rk, r]) for rk in ls])
            rxns = [ri(rk) if subT[rk, r] >= 0 else rev(ri(rk)) for rk in ls]
            dicto[subname(r, rxns)] = {'reacs': rxns,
                'coeffs': [abs(subT[rk, r]) for rk in ls],
                'full': [(ri(rk), subT[rk, r]) for rk in ls], 
                'bounds': (max(0, rc.lower_bound), max(0, rc.upper_bound))}
            dicto_keys.append(subname(r, rxns))
            if rc.reversibility:
                rxns = [ri(rk) if subT[rk, r] >= 0 else rev(ri(rk)) for rk in ls]
                dicto[subnamerev(r, rxns)] = {'reacs': rxns,
                    'coeffs': [abs(subT[rk, r]) for rk in ls],
                    'full': [(ri(rk), -subT[rk, r]) for rk in ls],
                    'bounds': (abs(min(0, rc.upper_bound)), abs(rc.lower_bound))}
            #improve: fva function flip should already handle the weird cases of negative only bounds         
    return dicto, dicto_keys
     

def make_metatool_data(df, model, dicto_keys):
    lines = ''
    nb_ext = 1
    for i, column in enumerate(df.columns):
        dfc = df[column]
        ftime=True
        lines += dicto_keys[i] + ' : '
        for k, v in dfc[dfc < 0].items():
            if ftime:
                ftime=False
            else:
                lines += ' + '
            lines += str(-v) + ' ' + k
        if ftime:
            lines += f'ext_{nb_ext}'
            nb_ext += 1
        lines += ' = ' if model.reactions.get_by_id(column).reversibility else ' => '
        ftime=True
        for k, v in dfc[dfc > 0].items():
            if ftime:
                ftime=False
            else:
                lines += ' + '
            lines += str(v) + ' ' + k
        if ftime:
            lines += f'ext_{nb_ext}'
            nb_ext += 1
        lines += '\n'
    lines = '-METEXT\n' + ' '.join(['ext_' + str(i+1) for i in range(nb_ext-1)]) + '\n' + '\n-CAT\n' + lines
    return lines
    

def compress_model(model, target_reactions):
    """
    Main work function
    Compresses the model
    Converts it into a MetaTool file
    """
    original_model = model.copy()
    compressed_model = model.copy()
        
    model = fva_model(original_model, compressed_model)
    constrained_model = original_model.copy()
    
    subT = efmtool4cobra.compress_model_sympy(model) # subT is a matrix for conversion of flux vectors between the full and compressed model
    reduced_matrix = cobra.util.array.create_stoichiometric_matrix(model, array_type='lil')
    rev_reduced_matrix = [int(r.reversibility) for r in model.reactions]

    print('-- Compressed network --')

    names = lambda l: list(map(lambda x: x.id, l))
    df = pandas.DataFrame.sparse.from_spmatrix(reduced_matrix, index=names(model.metabolites), columns=names(model.reactions))
    lines = get_reactions_equations(df, model)
    
    rxn_in_sub = [numpy.where(subT[:, i])[0] for i in range(subT.shape[1])]
            
    dicto, dicto_keys = make_compression_dict(model, original_model, subT, rxn_in_sub)
    
    objectives = []
    for rsub in dicto.keys():
        reacs = rsub.split('*')
        if set(reacs).intersection(set(target_reactions)):
            objectives.append(rsub)
    
    if objectives:
        target_reactions = objectives
    
    metatool = make_metatool_data(df, model, dicto_keys)
            
    return metatool, dicto, target_reactions



def compress_mnet(mnet: MetaNetwork, dict_file="", target_reactions=[]):
    """
    Compresses the MetaNetwork into another MetaNetwork
    Uses parsing functions for convenience
    Please check mparser parsing parameters, especially for SBMLs
    """
    with tempfile.NamedTemporaryFile('w', delete=False) as fp:
        format_sbml(mnet, fp.name)
        model = cobra.io.read_sbml_model(fp.name)
    metatool_model, subset_dict, target_reactions = compress_model(model, target_reactions)
    if dict_file:
        with open(dict_file, 'w') as f:
            pprint.pprint(subset_dict, stream=f, compact=True)
    with tempfile.NamedTemporaryFile('w', delete=False) as fp:
        fp.write(metatool_model); fp.close()  
        compressed_mnet = txt_read_network(fp.name)
    return compressed_mnet, target_reactions

