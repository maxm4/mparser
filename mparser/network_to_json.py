# -*- coding: utf-8 -*-

"""
network_to_json.py

Converts a metabolic network to a JSON file adapted for reading in Prolog

"""

from .meta_network import MetaNetwork
from .parsehelper import add_quotes 
import numpy as np
import json
    
def format_json(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to a JSON file usable notably in Prolog
    
    Params:
        network: MetaNetwork object
        out_fname: JSON format output file
    """
    
    json_dict = {'reactions': network.split_reactions,
                 'metabolites': network.metint,
                 'smatrix': network.rsmatrix.tolist(),
                 'rmatrix': network.rmatrix.tolist(),
                 'stoich': [list(tup) for tup in network.stoichiometry],
                 'neighbours': abs(np.sign(network.rsmatrix.T @ network.rsmatrix, dtype=int, casting='unsafe')).tolist()}#[(  for r in network.reactions]}

    with open(out_fname, 'w') as f:
        json.dump(json_dict, f)
    
