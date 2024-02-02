# -*- coding: utf-8 -*-

"""
network_to_bdd.py

Converts a metabolism network to a BDD

"""

from .meta_network import MetaNetwork
from .parsehelper import add_quotes, add_rev
import numpy as np
import pickle

def write_rrules(network, s, p, revs):
    print(f'{s} => {p}')
    if (s in revs) and (p in revs):
        print(f'{add_rev(p)} => {add_rev(s)}')

def format_bdd(network:MetaNetwork, out_fname):
    """
    Formats a metabolism network file to a BDD
    
    Params:
        network: MetaNetwork object
        out_fname: Hypergraph format output file
    """
    revs = network.reversibles
    print(network.idx2imet)
    print(network.idx2reac)

    rules = []
    svars = network.reactions
    for r in revs:
        svars.append(add_rev(r))
        rl = f'~ {r} | ~ {add_rev(r)}'
        rules.append(rl)
        print(rl)
    svars = sorted(svars)

    for j, r in enumerate(network.matrix):
        print(network.idx2imet[j], ':') 
        prods = np.where(r > 0)[0]
        prods = list(map(lambda x: network.idx2reac[x], prods))
        prods_rev = list(filter(lambda x: x in revs, prods))
        subst = np.where(r < 0)[0]
        subst = list(map(lambda x: network.idx2reac[x], subst))
        subst_rev = list(filter(lambda x: x in revs, subst))
        all_subst = subst + prods_rev
        for s in all_subst:
            # exclude s from all products and rev subs
            prods_ = list((set(prods) - set([s])))
            subst_rev_ = list(map(add_rev, (set(subst_rev) - set([s]))))
            djn = ' | '.join(prods_ + subst_rev_)
            head = s if s not in prods_rev else add_rev(s)
            if djn:
                rl = f'{head} => ({djn})'
                rules.append(rl)
                print(rl)

    svars = list(map(lambda x: x.replace('-', '_'), svars))
    rules = list(map(lambda x: x.replace('-', '_'), rules))
    rules = sorted(rules, key=lambda x: len(x), reverse=True)

    sv_fname = out_fname.replace('.dddmp', '.p')
    with open(sv_fname, 'wb') as f:
        pickle.dump(svars, f)

    try:
        from dd.cudd import BDD
        # Compute ROBDD using dd module
        compute_robdd = True # change to False to avoid computing it again
        if compute_robdd:
            bdd = BDD()
            bdd.declare(*svars)
            flist = []
            v = None
            for r in rules:
                print(r)
                u = bdd.add_expr(r)
                flist.append(u)
                if v is not None:
                    u &= v
                v = u
            bdd.dump(out_fname, [u])
    except ImportError:
        pass

    print(f'Created file {sv_fname}') # svars
    print(f'Created file {out_fname}') # bdd
    
