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


@click.command()
@click.option('--pk')
def main(pk):

    node = load_node(pk)
    descriptors = node.inputs.descriptor_names.get_list()
    species_definitions = node.inputs.species_definitions
    facet = species_definitions['s']['site_names'][0]
    surface = node.inputs.surface_names.get_list()[0]

    production_rate = node.outputs.production_rate_map.get_list()

    with open(f'output/node_{pk}_surface_{surface}_facet_{facet}.json',
              'w') as handle:
        json.dump(production_rate, handle)


if __name__ == '__main__':
    create_output_directory()

    main()
