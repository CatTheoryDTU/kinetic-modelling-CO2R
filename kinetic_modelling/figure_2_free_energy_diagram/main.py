# -*- coding: utf-8 -*-
# pylint: disable-all
"""
Script to plot Figure 2 of the manuscript.
More specifically, the output of this script has:
    1. Computational free energy diagram for Au, FeNC and NiNC
        - at standard conditions of temperature and pressure
        - with the potentials -0.6, -0.8, -1.0
        - with the pH = 2
    2. Replotting experimental data from Wen and Au
        - Stored in xls files
"""
import argparse
from glob import glob
import csv
import json
import string  # pylint: disable=syntax-error
from pprint import pprint
from pathlib import Path
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
import matplotlib.pyplot as plt
from experimental import plot_experimental_data
from computational_panel import FreeEnergyDiagram, plot_computational_diagram
from molecule import plot_molecule
from kinetic_modelling.plot_params import get_plot_params_Arial as get_plot_params  # pylint: disable=import-error, no-name-in-module


def cli_parse():
    """Parse the inputs from the command line, defaults are those used in the paper.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--database_folder',
                        default='../databases/',
                        help='Database folder')
    parser.add_argument('--potential', default=[-0.6, -0.8, -1.], nargs='*', type=float, \
                            help='SHE potential')
    parser.add_argument('--ph', default=2., type=float, help='pH')
    parser.add_argument('--referencedb_name',
                        default='reference_database/gas_phase.db')
    parser.add_argument('--ninc_experiment',
                        default='inputs/pH_effect_NiNC.xls')
    parser.add_argument('--fenc_experiment',
                        default='inputs/pH_effect_FeNC.xls')
    parser.add_argument('--fenc_xilehu_experiment',
                        default='inputs/pH_effect_FeNC_XileHu.xls')
    parser.add_argument('--gold_experiment',
                        default='inputs/pH_effect_Gold.xls')
    parser.add_argument('--copc_experiment',
                        default='inputs/pH_effect_CoPc.xls')
    parser.add_argument('--molecular_database',
                        default='inputs/molecule_CO2R.db')
    return parser.parse_args()


def main():  # pylint: disable=too-many-locals, too-many-statements
    """Main function that plots the free energy diagrams.
    """

    # Details of the figure
    fig = plt.figure(figsize=(20, 12.5))  # pylint: disable=invalid-name
    gs = fig.add_gridspec(8, 4)  # pylint: disable=invalid-name
    ax_Au = fig.add_subplot(gs[4:, 0])  # pylint: disable=invalid-name
    ax_Fe = fig.add_subplot(gs[4:, 1])  # pylint: disable=invalid-name
    ax_Ni = fig.add_subplot(gs[4:, 2])  # pylint: disable=invalid-name
    ax_Co = fig.add_subplot(gs[4:, 3:])  # pylint: disable=invalid-name
    cax_Au = fig.add_subplot(gs[2:4, 0])  # pylint: disable=invalid-name
    cax_Fe = fig.add_subplot(gs[2:4, 1])  # pylint: disable=invalid-name
    cax_Ni = fig.add_subplot(gs[2:4, 2])  # pylint: disable=invalid-name
    cax_Co = fig.add_subplot(gs[2:4, 3:])  # pylint: disable=invalid-name
    axf_Au = fig.add_subplot(gs[0:2, 0])  # pylint: disable=invalid-name
    axf_MNC = fig.add_subplot(gs[0:2, 1:3])  # pylint: disable=invalid-name
    axf_CoPc = fig.add_subplot(gs[0:2, 3])  # pylint: disable=invalid-name

    ax = [cax_Au, cax_Fe, cax_Ni, cax_Co, ax_Au, ax_Fe, ax_Ni, ax_Co]  # pylint: disable=invalid-name

    # Parse the required information from the command line
    parser = cli_parse()
    databases = glob(parser.database_folder + '/*.db')

    # Store relevant information in this dict
    data = {}

    # Experimental data for different experiments
    # Here fit_min refers to the minimum potentials
    # fit_max refers to the maximum potential at which
    # the Tafel slope was fit.
    # Refer to the main text for more details about what these
    # values mean
    plot_experimental_data(parser.gold_experiment,
                           ax_Au,
                           r'pH independent',
                           r'Au(pc)',
                           jmol_colors[atomic_numbers['Au']],
                           fit_min=-0.75,
                           fit_lim=-1.05)
    plot_experimental_data(parser.fenc_experiment,
                           ax_Fe,
                           r'pH independent',
                           r'FeNC',
                           jmol_colors[atomic_numbers['Fe']],
                           fit_min=-0.5,
                           fit_lim=-0.85)
    plot_experimental_data(parser.fenc_xilehu_experiment,
                           ax_Fe,
                           r'pH independent',
                           r'FeNC',
                           jmol_colors[atomic_numbers['Fe']],
                           fit_min=-0.5,
                           fit_lim=-0.85,
                           marker='v',
                           plot_line=False)
    ax_Fe.annotate(r'80 $\frac{\mathregular{mV}}{\mathregular{dec}}$',
                   xy=(0.75, 0.8),
                   xycoords='axes fraction')
    plot_experimental_data(parser.ninc_experiment,
                           ax_Ni,
                           r'pH dependent',
                           r'NiNC',
                           jmol_colors[atomic_numbers['Ni']],
                           fit_lim=-1.05,
                           fit_min=-0.75,
                           fit_all=False)
    plot_experimental_data(parser.copc_experiment,
                           ax_Co,
                           r'pH dependent',
                           r'CoPc',
                           jmol_colors[atomic_numbers['Co']],
                           fit_lim=-1.,
                           fit_min=-0.5,
                           fit_all=False,
                           pH_material='Co')

    ## Get the Free energy diagram
    for potential in parser.potential:
        method = FreeEnergyDiagram(dbnames=databases,\
                                   refdbname=parser.referencedb_name,\
                                   potential=potential,
                                   pH=parser.ph)
        method.main()

        data[potential] = method.diagram
        writeout = method.writeout
        writeout_zero = method.writeout_zero
        explicit_charge = method.explicit_charge
        E0 = method.E0  # pylint: disable=invalid-name

        filename = 'output/catmap_potential_%1.2f.txt' % potential
        with open(filename, 'w') as handle:
            csvwriter = csv.writer(handle, delimiter='\t')
            for row in writeout:
                csvwriter.writerow(row)

    # do the plot for the molecular part
    # this conforms with the new way of doing finite difference
    # using the newer implementation, which is why it is a new class
    plot_molecule([parser.potential[1]], parser.ph, parser.molecular_database,
                  cax_Co, method.references, method.references_E)

    # save the catmap input file
    filename = 'output/catmap_potential_pzc.txt'
    with open(filename, 'w') as handle:
        csvwriter = csv.writer(handle, delimiter='\t')
        for row in writeout_zero:
            csvwriter.writerow(row)

    print('------------')
    print('Charge from the finite difference method')
    pprint(explicit_charge)

    with open('../databases/explicit_charge.json', 'w') as handle:
        json.dump(explicit_charge, handle)

    with open('../databases/zero_charge_energies.json', 'w') as handle:
        json.dump(E0, handle)

    plot_computational_diagram(data, [cax_Au, cax_Fe, cax_Ni],
                               SAC_potential=-0.8)

    ## Label the diagram
    alphabet = list(string.ascii_lowercase)
    for i, ax_ in enumerate(ax):
        ax_.annotate(alphabet[i] + ')',
                     xy=(0.05, 0.87),
                     xycoords='axes fraction',
                     fontsize=20)

    # Images to be added to the figure
    arr_image = plt.imread('input_images/Au27.020.png', format='png')
    axf_Au.imshow(arr_image)
    axf_Au.axis('off')

    arr_image = plt.imread('input_images/defects.png', format='png')
    axf_MNC.imshow(arr_image)
    axf_MNC.axis('off')

    arr_image = plt.imread('input_images/CoPc.png', format='png')
    axf_CoPc.imshow(arr_image)
    axf_CoPc.axis('off')

    fig.tight_layout()
    fig.savefig('output/figure1.pdf')


if __name__ == '__main__':
    Path('output').mkdir(exist_ok=True)
    Path('output_si').mkdir(exist_ok=True)
    get_plot_params()
    main()
