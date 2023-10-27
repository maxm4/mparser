# -*- coding: utf-8 -*-

"""
__init__.py

Parser module
Reads a metabolic network
Converts a network to ASP

"""

import sys

if sys.version_info[:2] > (3,5):
    from .conversion import InputFormatType, OutputFormatType, convert, read_network
    from .network_to_sbml import format_sbml
    from .network_to_efmtool import format_efmtool
    from .network_to_asp import format_asp, write_asp, append_asp
    from .network_to_metatool import format_metatool
    from .meta_network import MetaNetwork

from . import parsehelper

__all__ = ['MetaNetwork', 'read_network', 'parsehelper', 'format_sbml',
           'format_asp', 'InputFormatType', 'OutputFormatType', 'write_asp', 
           'append_asp', 'format_efmtool', 'format_metatool']

