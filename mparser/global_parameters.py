"""
reading_parameters.py

Parameters for reading SBML files and other formats
Include the definition of external metabolites

"""

### As of now, this includes SBML Reading Only ###
# Lambda functions used to check if a species is external
# To be used accordingly depending of the SBML Version and Format
boundary_condition = lambda species: species.getBoundaryCondition() # very consistent method, toy biofuel
external_compartment = lambda species: species.getCompartment() == ext_cp_name
boundary_name = lambda species: species.getId().endswith(ext_suffix) # try this with orth_core ?
ext_suffix = '_b' # or '_e' or '_p', note: in iAF1260 '_e' still is not the final external meta delimiter
ext_cp_name = "e" # or "C_e", etc., this is another parameter
remove_null_reactions = False # if True, completely excludes null reactions from MetaNetwork object
metext_func = boundary_name # external_compartment # condition # boundary_name
# boundary_name with ext_suffix: '_b' is the recommended method with all Cobra models
# No external metabolites should be detected, and that's normal
# Indeed, in standard Cobra representation, all Cobra metabolites are assumed internal
# Including metabolites from the external compartments
### End of values used by sbml_to_network.py ###
