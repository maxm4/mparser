# mparser

Module for converting metabolic networks into Answer Set Programming (ASP) and more

Pip requirements are in `requirements.txt`.

Companion module of *[aspefm](https://github.com/maxm4/aspefm)*.

## Main scripts :

`mparser_cli.py`: Allows conversion of networks between various formats, command line interface

`mparser_gui.py`: Allows conversion of networks between various formats, graphical user interface (requires PySimpleGui)

## Additional scripts :

`aspsmatrix.py`: Extracts and pickles the matrix from the mparser-read matrix, also contains many legacy functions

`parse_and_check_output.py`: More advanced version of `aspefm/parse_output.py` using the kernel test on the aspsmatrix matrix to check validity of solutions

`thermo_csv.py`: Short script converting a csv of keq values into an ASP file. See `data` folder for thermo data.

`prepare_efm_checker.py`: Script preparing the input file for `aspefm/extensions/efm_checker.py` based on the aspsmatrix matrix.

`add_sbml_boundaries.py`: Script converting COBRA SBML files into non-COBRA SBML files with added boundary metabolites in the most outer compartment.

Some of these scripts contain legacy functions, which may not be still functional today.

## The module :

Implements conversion from METATOOL format and SBML format to ASP and more. Found in the `mparser` folder.

The file `global_parameters.py` is meant as a file gathering all global parameters that are not given in argument to `mparser_cli.py` and `mparser_gui.py`.

SBML to network reading parameters are by default adapted to COBRA files, which are not necessarily representative of every possible SBML.

## Example usages :

```python mparser_cli.py METATOOL data/ppp.txt ASP data/ppp.lp4 ```

```python mparser_cli.py METATOOL data/ppp.txt SBML data/ppp.xml ```

```python mparser_cli.py SBML data/ppp.xml ASP data/ppp_dual.lp4 --to-dual-mcs --target-reactions R_ALD```

```python prepare_efm_checker.py $cobrafile.xml $cobrafile.efmc```

```python parse_and_check_output.py $file.txt $file.efms --check $sbmlNetwork --pickle $file.pkl --sbml```
