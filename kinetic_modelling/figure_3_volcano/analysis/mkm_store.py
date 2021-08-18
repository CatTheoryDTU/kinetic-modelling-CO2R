# -*- coding: utf-8 -*-
# pylint: disable-all
"""

For a given CalcJobNode store the production map

Needs the following inputs

* Node of the CalcJobNode where CatMAP was the parser

"""

import sys
import os
import numpy as np
from pathlib import Path
from useful_classes import bcolors
import matplotlib
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from pprint import pprint
from useful_functions import get_fit_from_points
import click
from ase.thermochemistry import HarmonicThermo
from matplotlib import ticker
from useful_functions import create_output_directory
from plot_params import get_plot_params
import json
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import argparse


def cli_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pks',
                        type=int,
                        nargs='*',
                        help='Any number of input pk')
    return parser.parse_args()


def main():

    parser = cli_parser()

    pks = parser.pks
    data_tot = {}

    for pk in pks:

        node = load_node(pk)
        descriptors = node.inputs.descriptor_names.get_list()
        species_definitions = node.inputs.species_definitions
        facet = species_definitions['s']['site_names']
        surfaces = node.inputs.surface_names.get_list()
        potential = node.inputs.voltage.value
        pH = node.inputs.pH.value
        coverage_map = node.outputs.coverage_map.get_list()
        production_rate = node.outputs.production_rate_map.get_list()
        energy_file = node.inputs.energies.get_content()

        data = {}
        data['facet'] = facet
        data['surfaces'] = surfaces
        data['descriptors'] = descriptors
        data['potential'] = potential
        data['pH'] = pH
        data['coverage_map'] = coverage_map
        data['production_rate'] = production_rate
        data['energy_file'] = energy_file
        data_tot[pk] = data

    with open('aiida_output/kinetic_model_data.json', 'w') as handle:
        json.dump(data_tot, handle, indent=4)


if __name__ == '__main__':
    create_output_directory('aiida_output')
    main()
