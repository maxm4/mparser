# -*- coding: utf-8 -*-

"""
prepare_efm_checker.py

Module computing an EFMChecker config file
From a SBML metabolic network or AspSmatrix in Pickle format
Using module aspsmatrix

"""

from types import SimpleNamespace
from aspsmatrix import *
import pickle
import ast

def compute_efm_checker(infile, outfile, litfile, pkl):
    """
    Computes an EFMChecker config file
    From a SBML metabolic network
    Using module aspsmatrix

    """
    if not pkl:
        collec = sbml_file(infile)
        aspsm = AspSmatrix(collec)
    else:
        with open(infile, 'rb') as f:
            metanw = pickle.load(f)
        aspsm = AspSmatrix(metanw, from_files=False)

    rindex = aspsm.reactions_index()
    matrix = aspsm.matrix
    neighb = compute_neighbours(aspsm)
    literals = {}

    if litfile:
        with open(litfile, 'r') as f:
            literals = ast.literal_eval(f.readline().strip())

    config = SimpleNamespace(rindex=rindex, matrix=matrix, neighbours=neighb, literals=literals)

    with open(outfile, 'wb') as f:
        pickle.dump(config, f)

if __name__== "__main__":
    parser = ArgumentParser()
    parser.add_argument('infile', metavar='input.file', help='Input file name')
    parser.add_argument('outfile', metavar='output.file', help='Output file name')
    parser.add_argument('--litfile', metavar='literals.file', help='File containing literals')
    parser.add_argument('--pickle', action='store_true', help='If network file is in Pickle format')
    opts = parser.parse_args()
    compute_efm_checker(opts.infile, opts.outfile, opts.litfile, opts.pickle)


