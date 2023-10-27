# -*- coding: utf-8 -*-
"""
dualize_mcs.py

Converts MetaNetwork object to dual

"""

from copy import deepcopy
import numpy as np
from .meta_network import MetaNetwork

def dual_mcs(mnet, target_reactions=[], irr_reactions=False):
    """
    Converts the network to a dual network for computing
    Minimal Cut Sets
    
    Params:
        target_reactions: target reaction names
        irr_reactions: if irreversible mock reactions should be included,
                       following the original Ballerstein et al formulation
    """
    # Indicating conversion to dual network
    mnet.old = deepcopy(mnet)
    mnet.dual_mcs = True
    # Stoichiometry matrix shenanigans
    mnet.matrix = mnet.matrix.T
    reacs = np.eye(mnet.nb_reactions)
    if irr_reactions:
        irrevs = -reacs.copy()
        idx_rev = [mnet.reac2idx[r] for r in mnet.reversibles] 
        irrevs = np.delete(irrevs, idx_rev, axis=1)
    targets = -np.ones((mnet.nb_reactions, 1))
    untargeted = mnet._reactions_s ^ set(target_reactions)
    idx_ut = [mnet.reac2idx[r] for r in untargeted] 
    targets[idx_ut] = 0
    if irr_reactions:
        mnet.matrix = np.concatenate((mnet.matrix, reacs, irrevs, targets), axis=1)
    else:
        mnet.matrix = np.concatenate((mnet.matrix, reacs, targets), axis=1)
    # Dual reaction names definition
    namecol = lambda x: 'mcs_' + x
    addtarg = lambda x: x + '_tgt'
    mnames = list(map(namecol, mnet.metint))
    rnames = list(map(namecol, mnet.reactions))
    xnames = list(map(namecol, mnet.transporters))
    if irr_reactions:        
        addirr = lambda x: x + '_irr'
        irnames = list(map(namecol, np.delete(mnet.reactions, idx_rev, axis=0)))
        irnames = list(map(addirr, irnames))
    tnames = namecol('_'.join(target_reactions))
    tnames = [addtarg(tnames)]
    # Transfer references to work with 
    metext = []
    metabolites = mnet.reactions
    if irr_reactions:
        reactions = mnames + rnames + irnames + tnames
        reversibles = mnames + rnames
    else:
        fwdrnames = sorted(list(set(rnames) - set(list(map(namecol, mnet.irreversibles))) ))
        reactions = mnames + rnames + tnames
        reversibles = mnames + fwdrnames
    stoichiometry = []
    # Reactions of interest
    mnet.extras = {} # new extra field grouping -- name must be singular
    mnet.extras['interest'] = rnames
    mnet.extras['stoichreac'] = mnames
    mnet.extras['mirrev'] = sorted(list(mnet.irreversibles))
    if irr_reactions:
        mnet.extras['parallel'] = irnames
    mnet.transporters = xnames
    # Additional dicts to help with direct association
    mnet.idx2imet = {i: m for i, m in enumerate(metabolites)} 
    mnet.idx2reac = {i: r for i, r in enumerate(reactions)}
    # Stoichiometry recreation
    for i, row in enumerate(mnet.matrix):
        for j, case in enumerate(row):
            if case != 0:
                stoichiometry.append((mnet.idx2imet[i], mnet.idx2reac[j], case))
    # Finalize conversion process by recreating object
    mnet.__init__(metext, metabolites, reactions, stoichiometry, reversibles)
    return mnet

