# -*- coding: utf-8 -*-

"""
network_to_hypergraph.py

Converts a metabolism network to an hypergraph in PaToH format

"""

from .meta_network import MetaNetwork
from .parsehelper import add_quotes 
import numpy as np
    
def format_hypergraph(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to a standard hypergraph format
    
    Params:
        network: MetaNetwork object
        out_fname: Hypergraph format output file
    """
    hypergraph = ""

    netcount = 0
    print(network.idx2imet)
    print(network.idx2reac)
    for r in network.matrix.T:
        indices = np.nonzero(r)[0]
        netcount += len(indices)
        indices = " ".join(list(map(str,indices)))
        hypergraph += f"{indices}\n"

    firstline = f"0 {network.nb_metabolites} {network.nb_reactions} {netcount}\n"
    hypergraph = firstline + hypergraph

    f_out = open(out_fname, "w")
    f_out.write(hypergraph)
    f_out.close()  
    
