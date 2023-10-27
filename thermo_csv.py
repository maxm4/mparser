# -*- coding: utf-8 -*-

"""
thermo_csv.py

Specific CSV parsing for Keq equilibrium constants
Parses the first column as reaction names, the third column as log Keq
Writes ASP constraints with the corresponding values

"""

import sys, csv    

def thermo_csv_script():
    """
    
    Script execution when this file is run as standalone
    Parses the first column as reaction names, the third column as log Keq
    Writes ASP constraints with the corresponding values
    
    """
    if len(sys.argv) > 2:
        in_name = sys.argv[1]
        out_name = sys.argv[2]
    else:
        print("Please enter file names - ")
        in_name = input("CSV file: ")
        out_name = input("ASP file: ")
        
    asp = "" 
    withquotes = lambda w: f'"{w}"'
    withdot = lambda w: w.replace(',','.')
    with open(in_name, "r", newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            assert(len(row) > 2)
            name = withquotes(row[0])
            name_rev = withquotes(row[0] + "_rev")
            if ',' in row[2]:
                value = withquotes(withdot(row[2]))
                neg_value = withquotes(str(-float(withdot(row[2]))))
            else: # only adapted to french xls locale...
                value = str(int(row[2]))
                neg_value = str(-int(row[2]))
            asp += "keq(" + name + "," + value + ").\n"
            asp += "keq(" + name_rev + "," + neg_value + ").\n"
            
    with open(out_name, "w") as g:
        g.write(asp)
    

def read_csv_values(fname):
    """
    Reads a csv file containing names and values
    Parses the first column as names and the second column as values
    
    Returns:
        names: list of names
        values: list of values
    """
    names = []
    values = []
    with open(fname, "r", newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            assert(len(row) > 1)
            if ',' in row[1]:
                row[1] = row[1].replace(',','.') # only adapted to french xls locale...
            names.append(row[0])
            values.append(float(row[1]))
    return names, values


if __name__== "__main__":
    thermo_csv_script()

