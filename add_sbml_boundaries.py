from argparse import ArgumentParser
import cobra
import getpass
from datetime import datetime

ext_suffix = '_e' # or '_b' or '_p', note: in iAF1260 '_e' still is not the final external meta delimiter
new_ext_suffix = '_b' # or '_xt', new suffix to distinguish from the previous suffix (ie. ext_suffix)
replace_suffix = lambda n: n[:-len(ext_suffix)] + new_ext_suffix
add_suffix = lambda n: n + new_ext_suffix
add_history = False
boundary_correct = True

def correct_boundary(rb, apply_suffix=replace_suffix):
    """
    Corrects a certain type of COBRA Model boundary reaction
    Adds a boundary metabolite to the COBRA Model
    
    Params:
        rb: string representing the boundary reaction
        apply_suffix: lambda function determining what should be done with the suffix
        
    Returns:
        rb: string representing the updated boundary reaction
        
    Raises:
        AssertionError: if the COBRA Model boundary reaction is badly formed
    """
    tab = rb.split()
    assert(len(tab) == 2)
    assert(tab[1] == '-->' or tab[1] == '<=>' or tab[1] == '<--')
    return tab[0]+ ' ' + tab[1] + ' ' + apply_suffix(tab[0])



def cobra_to_sbml(infname, outfname):
    """
    Converts a COBRA SBML to an up-to-date SBML with standard boundary metabolites
    Uses Python library COBRA
    
    Params:
        infname: name of the input COBRA SBML file
        outfname: name of the output SBML file
    """
    m = cobra.io.read_sbml_model(infname)
    if boundary_correct:
        for x in m.boundary:
            x.reaction = correct_boundary(x.reaction)
        no_comp = [x for x in m.metabolites if x.compartment is None]
        for x in no_comp:
            x.compartment = 'b'
    cobra.io.sbml.write_sbml_model(m, outfname)
    
    if add_history:
        username = getpass.getuser()
        today = datetime.today().strftime("%Y-%m-%d at %H:%M:%S")
        with open(outfname, "r+") as fp:
            lines = fp.readlines()
            lines.insert(0, f"<!-- last edited by {username} on {today} -->\n")
            fp.seek(0)
            fp.writelines(lines)

if __name__== "__main__":
    parser = ArgumentParser()
    parser.add_argument('infile', help='Input filename: Original COBRA SBML')
    parser.add_argument('outfile', help='Output filename: Updated SBML')
    opts = parser.parse_args()
    cobra_to_sbml(opts.infile, opts.outfile)
