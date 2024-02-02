# -*- coding: utf-8 -*-

"""
conversion.py

Module containing format types specifications and linking the different conversion methods together
Converts a file specificating a metabolic network from a format to another

"""

from enum import Enum

from .asp_to_network import asp_read_network
from .dualize_mcs import dual_mcs
from .compress_model import compress_mnet
from .meta_network import MetaNetwork
from .network_to_asp import format_asp
from .network_to_sbml import format_sbml
from .network_to_efmtool import format_efmtool
from .network_to_metatool import format_metatool
from .network_to_scrumpy import format_scrumpy
from .network_to_hypergraph import format_hypergraph
from .network_to_bdd import format_bdd
from .network_to_json import format_json
from .scrumpy_to_network import scrumpy_read_network
from .metatool_to_network import txt_read_network
from .sbml_to_network import sbml_read_network
import pickle

class InputFormatType(Enum):  
    """ 
    
    InputFormatType enum 
    Defines constants corresponding to allowed input format types
    
    """
    TXT = 'METATOOL'
    SBML = 'SBML'
    SPY = 'SCRUMPY'
    ASP = 'ASP'
    PKL = 'PICKLE'

    def __str__(self):
        return self.value


class OutputFormatType(Enum):  
    """ 
    
    OutputFormatType enum 
    Defines constants corresponding to allowed output format types
    
    """
    TXT = 'METATOOL'
    SBML = 'SBML'
    ASP = 'ASP'
    EFM = 'EFMTOOL'
    PKL = 'PICKLE'
    SPY = 'SCRUMPY'
    HYP = 'GRAPH'
    BDD = 'BDD'
    JSON = 'JSON'

    def __str__(self):
        return self.value
    
    
def pickle_read_network(in_fname):
    """
    Reads a metabolic network with pickle
    
    Params:
        in_fname: input metabolism network file name
    
    """    
    with open(in_fname, 'rb') as infile:
        return pickle.load(infile)


def read_network(in_fname, in_format):
    """
    Reads a metabolism network from a file depending of the given format
    Returns a MetaNetwork object
    
    Params:
        in_fname: input metabolism network file name
        in_format: input metabolism network format
        
    Returns:
        network: MetaNetwork object
    """
    if in_format == InputFormatType.TXT:
        return txt_read_network(in_fname)
    elif in_format == InputFormatType.SBML:
        return sbml_read_network(in_fname)
    elif in_format == InputFormatType.SPY:
        return scrumpy_read_network(in_fname)
    elif in_format == InputFormatType.ASP:
        return asp_read_network(in_fname)
    elif in_format == InputFormatType.PKL:
        return pickle_read_network(in_fname)
    else:
        raise NotImplementedError # unsupported format

        
def format_pickle(network:MetaNetwork, out_name):
    """
    Dumps a metabolic network with pickle
    
    Parameters:
        network: MetaNetwork instance
        out_name: output file name
    
    """    
    with open(out_name, 'wb') as outfile:
        pickle.dump(network, outfile)

    
def convert(in_format:InputFormatType, in_name, 
            out_format:OutputFormatType, out_name,
            compression:bool, dict_file_name:str,
            to_dual_mcs, target_reactions, ballerstein):
    """ 
    
    Converts the metabolic network in the input file from a format to another
    Outputs the result in the output file
    
    Parameters:
        in_format: InputFormatType instance
        in_name: input file name
        out_format: OutputFormatType instance
        out_name: output file name
        to_dual_mcs: converts to a dual network for computing Minimal Cut Sets
        ballerstein: Boolean indicating if Ballerstein formulation for MCSs
    """
    meta_network = read_network(in_name, in_format)
    if compression:
        meta_network, target_reactions = compress_mnet(meta_network, dict_file_name, target_reactions)
    if to_dual_mcs:
        meta_network = dual_mcs(meta_network, target_reactions=target_reactions, irr_reactions=ballerstein)
    if out_format == OutputFormatType.ASP:
        format_asp(meta_network, out_name)
    elif out_format == OutputFormatType.SBML:
        format_sbml(meta_network, out_name)
    elif out_format == OutputFormatType.EFM:
        format_efmtool(meta_network, out_name)
    elif out_format == OutputFormatType.PKL:
        format_pickle(meta_network, out_name)
    elif out_format == OutputFormatType.TXT:
        format_metatool(meta_network, out_name)
    elif out_format == OutputFormatType.SPY:
        format_scrumpy(meta_network, out_name)
    elif out_format == OutputFormatType.HYP:
        format_hypergraph(meta_network, out_name)
    elif out_format == OutputFormatType.BDD:
        format_bdd(meta_network, out_name)
    elif out_format == OutputFormatType.JSON:
        format_json(meta_network, out_name)
    else:
        raise NotImplementedError # unsupported format
    
    
