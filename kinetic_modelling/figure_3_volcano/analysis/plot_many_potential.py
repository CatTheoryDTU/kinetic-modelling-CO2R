# -*- coding: utf-8 -*-
# pylint: disable-all
import click
from pprint import pprint
import json
from useful_functions import create_output_directory
import numpy as np
from ase.data import atomic_numbers
from pprint import pprint
from useful_functions import get_fit_from_points
from scipy.interpolate import griddata, interp2d, Rbf
import matplotlib.pyplot as plt
from ase.data.colors import jmol_colors
from matplotlib import ticker
from plot_params import get_plot_params
import xlrd
import string
from plot_kinetics_figure import get_electronic_energy_points, plot_map


def main():

    pk_potential = {
        -0.2: 12425,
        -0.4: 12454,
        -0.6: 12483,
        -0.8: 12384,
        -1.0: 12512,
        -1.2: 12541,
        # -1.4: 12570,
    }

    fig, axf = plt.subplots(3, 2, constrained_layout=True, figsize=(15, 15))

    kfiles = 'aiida_output/kinetic_model_data.json'
    with open(kfiles, 'r') as handle:
        data_tot = json.load(handle)

    ax = axf.ravel()

    for i, (potential, pk) in enumerate(pk_potential.items()):
        data = data_tot[str(pk)]
        data_points = get_electronic_energy_points(data['energy_file'],
                                                   data['surfaces'],
                                                   data['facet'][0])
        # if i == 0: plot_legend = True ; else: plot_legend = False
        plot_legend = True if i == 0 else False
        plot_map(fig=fig,
                 ax=ax[i],
                 maps=data['production_rate'],
                 descriptors=data['descriptors'],
                 points=data_points,
                 potential=data['potential'],
                 pH=data['pH'],
                 plot_single_atom=True,
                 cmapname='coolwarm',
                 plot_metal=True,
                 annotate_rate_limiting=False,
                 plot_legend=plot_legend)

    fig.savefig('output/potential_dependent_production.pdf')


if __name__ == '__main__':
    get_plot_params()
    main()
