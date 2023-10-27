# -*- coding: utf-8 -*-

"""
mparser_cli.py

Command Line Interface
Reads a metabolic network
Performs the requested conversion

"""

from mparser import OutputFormatType, InputFormatType
from mparser import convert
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument('informat', type=InputFormatType, choices=list(InputFormatType), help='Input format type') 
parser.add_argument('infile', help='Input file name')
parser.add_argument('outformat', type=OutputFormatType, choices=list(OutputFormatType), help='Output format type')
parser.add_argument('outfile', help='Output file name')
parser.add_argument('--target-reactions', help='Target reactions for computing Minimal Cut Sets', default=[])
parser.add_argument('--ballerstein', help='Enables Ballerstein formulation for Minimal Cut Sets', action='store_true')
parser.add_argument('--to-dual-mcs', action='store_true', help='Convert to dual network for computing Minimal Cut Sets')
opts = parser.parse_args()
in_format = opts.informat
in_name = opts.infile
out_format = opts.outformat 
out_name = opts.outfile
to_dual_mcs = opts.to_dual_mcs
ballerstein = opts.ballerstein
target_reactions = opts.target_reactions
if target_reactions: target_reactions = target_reactions.split(',')
convert(in_format, in_name, out_format, out_name, to_dual_mcs, target_reactions, ballerstein)
