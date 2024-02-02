# -*- coding: utf-8 -*-

"""
mparser_gui.py

Graphical User Interface
Reads a metabolic network
Performs the requested conversion

"""

from mparser import OutputFormatType, InputFormatType
from mparser import convert
import PySimpleGUI as sg
import sys

sg.theme('SystemDefaultForReal')
inputformats = list(InputFormatType)
outputformats = list(OutputFormatType)

layout = [[sg.Text('Input format', size=(15, 1)), sg.InputCombo(inputformats, inputformats[0])],
          [sg.Text('Input file', size=(15, 1)), sg.Input(), sg.FileBrowse()],
          [sg.Text('Output format', size=(15, 1)), sg.InputCombo(outputformats, outputformats[0])],
          [sg.Text('Output file', size=(15, 1)), sg.InputText(), sg.FileSaveAs()],
          [sg.Checkbox('Compress network (Enzyme subsets)', default=False)],
          [sg.Checkbox('To dual network (Minimal Cut Sets)', default=False)],
          [sg.OK(), sg.Cancel()] ]

window = sg.Window('Mparser', layout)
event, values = window.read()
window.close()

in_format = InputFormatType(values[0])
in_name = values[1]
out_format = OutputFormatType(values[2])
out_name = values[3]
to_dual_mcs = values[5]
compress = values[4]
ballerstein = False
target_reactions = []
compression_dict = str()

if compress: 
    layout = [[sg.Text('Compression dict', size=(15, 1)), sg.InputText(), sg.FileSaveAs()],
              [sg.OK(), sg.Cancel()]]
    window = sg.Window('Mparser', layout)
    event, values = window.read()
    window.close()
    compression_dict = values[0]
    

if to_dual_mcs:
    layout = [[sg.Text('Target reactions', size=(15, 1)), sg.InputText()],
              [sg.Text('* reactions should be separated by commas (,)', size=(50, 1))],
              [sg.Checkbox('Ballerstein dual formulation instead of von Kamp\'s', default=False)],
              [sg.OK(), sg.Cancel()]]

    window = sg.Window('Mparser', layout)
    event, values = window.read()
    window.close()
    target_reactions = values[0].split(',')
    ballerstein = values[1]
    
convert(in_format, in_name, out_format, out_name, compress, compression_dict, to_dual_mcs, target_reactions, ballerstein)
