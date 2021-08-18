# -*- coding: utf-8 -*-
# pylint: disable-all
import numpy as np
# from scipy.signal import hilbert
from scipy import signal
import matplotlib.pyplot as plt
from useful_functions import create_output_directory
from pprint import pprint
import numpy as np
from plot_params import get_plot_params
import click
import json
from scipy.signal import hilbert
from ase.db import connect
from useful_functions import get_fit_from_points
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
from scipy.integrate import quad, simps
from scipy.optimize import curve_fit
from ase.io import read
from ase.data import covalent_radii as radii
import string
from matplotlib.patches import Circle
import matplotlib.image as mpimg


def parsedb(results, database, sac=False, dbconfig={}):

    for row in database.select(*dbconfig):
        if sac:
            facet = row.vacancy_number.replace(
                'vacancy_', '') + '_' + row.dopant_number.replace(
                    'dopant_', '')
            metal = row.metal_dopant.replace('metal_dopant_',
                                             '').replace('_nonorth', '')
        else:
            facet = row.facets.replace('facet_', '')
            metal = row.sampling.replace('sampling_', '')
        state = row.states.replace('state_', '')
        results.setdefault(facet, {}).setdefault(metal, {}).setdefault(
            state, {})['dipole'] = row.dipole_field


def main():

    tmdbname = '../databases/transition_metal_vacuum.db'
    sacdbname = '../databases/single_atom_vacuum.db'

    results = {}
    parsedb(
        results,
        connect(tmdbname),
    )
    parsedb(results, connect(sacdbname), sac=True)

    # fig, ax = plt.subplots(1, 1, figsize=(9,4.5))
    fig = plt.figure(constrained_layout=True, figsize=(15, 8))
    gs = fig.add_gridspec(4, 7)
    ax = fig.add_subplot(gs[0:2, :4])
    axc = fig.add_subplot(gs[0:2, 4:])
    # axd = fig.add_subplot(gs[0:2,6:])

    axp = []
    for i in range(7):
        axp.append(fig.add_subplot(gs[2:, i]))
    for facet in results:
        if facet == '111':
            continue
        if '_' in facet:
            facetname = '(' + facet.replace('_', ',') + ')'
            color = 'tab:red'
        else:
            facetname = '(' + facet + ')'
            color = 'tab:blue'
        for metal in results[facet]:
            if metal == 'gas_molecules':
                continue
            elif metal == 'Al':
                continue

            xlabel = metal + facetname
            ax.plot(xlabel,
                    results[facet][metal]['CO']['dipole'],
                    'o',
                    color=color)
            ax.plot(xlabel,
                    results[facet][metal]['COOH']['dipole'],
                    'v',
                    color=color)
            ax.plot(xlabel,
                    results[facet][metal]['CO2_dos']['dipole'],
                    '*',
                    color=color)
    ax.tick_params(axis='x', labelrotation=90)
    ax.set_ylabel(r'$\mu$ / e$\AA$')
    ax.plot([], [], '*', color='k', label=r'CO$_2$*')
    ax.plot([], [], 'v', color='k', label='COOH*')
    ax.plot([], [], 'o', color='k', label='CO*')
    # ax.legend(loc='best', frameon=False)
    # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
    # mode="expand", borderaxespad=0, ncol=3, frameon=False, fontsize=14)
    # ax.legend( bbox_to_anchor=(1.04,1), borderaxespad=0, fontsize=12)
    ax.annotate('Transition Metals',
                xy=(0.2, 0.9),
                color='tab:blue',
                xycoords='axes fraction',
                fontsize=16)
    ax.annotate('MNC',
                xy=(0.7, 0.9),
                color='tab:red',
                xycoords='axes fraction',
                fontsize=16)
    # ax.legend(bbox_to_anchor=(1.04,0), loc="lower left", borderaxespad=0, fontsize=14, frameon=False)
    ax.legend(loc='best', frameon=False, fontsize=12)

    tm_files = '../databases/TM_dos.json'
    sac_files = '../databases/SAC_dos.json'
    gas_files = '../databases/Gas_dos.json'

    tm_data = json.load(open(tm_files, 'r'))
    sac_data = json.load(open(sac_files, 'r'))
    gas_data = json.load(open(gas_files, 'r'))

    arrow = {
        'Ag': [0.8, -2.2, 0, 0.6],
        'Au': [0.8, -2., 0, 1.],
        'Cu': [0.8, -2.2, 0, 1.2],
        'Pt': [0.8, -2.3, 0, 2.5],
        'Pd': [0.8, -2.5, 0, 3.0],
    }

    order_plots = {'Ag': 2, 'Au': 3, 'Cu': 4, 'Pd': 5, 'Pt': 6}

    for metal in tm_data:
        if metal == 'Al':
            continue
        for j, facet in enumerate(tm_data[metal]):
            if facet == '100':
                continue

            pdos_slab = np.array(
                tm_data[metal][facet]['pdos'][f'{metal}(d)']['+'][0])
            pdos_ads = np.array(
                tm_data[metal][facet]['pdos']['CO2(sp)']['+'][0])
            energies_ads = np.array(tm_data[metal][facet]['energies_ads'])
            energies_slab = np.array(tm_data[metal][facet]['energies_slab'])

            filled_indices = [
                i for i in range(len(energies_ads)) if energies_ads[i] <= 0.0
            ]
            filling = np.trapz(pdos_ads[filled_indices],
                               energies_ads[filled_indices]) / np.trapz(
                                   pdos_ads, energies_ads)
            print(metal, facet, filling)
            i = order_plots[metal]
            # axp[i].fill_between(pdos_slab,\
            #         energies_slab,\
            #         color='tab:purple', alpha=0.25)
            axp[i].plot(pdos_slab,\
                    energies_slab,\
                    color='tab:purple', alpha=0.3)
            axp[i].plot(pdos_ads,\
                    energies_ads,\
                    color='tab:green')
            axp[i].fill_between(pdos_ads[filled_indices],\
                    energies_ads[filled_indices],\
                    color='tab:green', alpha=0.25)
            axp[i].set_ylim([-5, 4])
            axp[i].annotate(
                r'$\mathregular{%s}\left (%s \right ) - \left ( \mathregular{d} \right )$'
                % (metal, facet),
                xy=(0.2, 0.01),
                xycoords='axes fraction',
                fontsize=14,
                color='tab:purple')
            axp[i].set_xticks([])
            axp[i].arrow(*arrow[metal],
                         width=0.05,
                         head_length=0.0,
                         head_width=0.0,
                         color='tab:green')
            axp[i].set_xlim([0.0, 1.2])
            if i != 0:
                axp[i].set_yticks([])
            else:
                axp[i].set_ylabel(r'$\epsilon - \epsilon_{f}$ / eV')
                # axp[i].plot([],[], color='tab:blue', lw=10, label=r'CO$_{2}*$')
                axp[i].annotate(r'$\mathregular{CO}_{2}^{*} \left ( \mathregular{s,p} \right)$', xy=(0.2,0.9),\
                     color='tab:green', xycoords='axes fraction', fontsize=14)
                axp[i].legend(loc='best', frameon=False, fontsize=14)
    i = 0
    for metal in sac_data:
        for value in sac_data[metal]:
            vacancy, dopant = value.split('_')
            if vacancy == '2' and dopant == '4':
                pass
            else:
                continue
            pdos_slab = []
            pdos_ads = []

            for spin in sac_data[metal][value]['pdos'][f'{metal}(d)']:
                pdos_slab.append(
                    sac_data[metal][value]['pdos'][f'{metal}(d)'][spin][0])
                pdos_ads.append(
                    sac_data[metal][value]['pdos']['CO2(sp)'][spin][0])

            energies_ads = np.array(sac_data[metal][value]['energies_ads'])
            energies_slab = np.array(sac_data[metal][value]['energies_slab'])

            summed_dos_slab = np.array(pdos_slab).sum(axis=0)
            summed_dos_ads = np.array(pdos_ads).sum(axis=0)
            axp[i].plot(summed_dos_ads, energies_ads, color='tab:green')
            filled_indices = [
                i for i in range(len(energies_ads)) if energies_ads[i] <= 0.0
            ]
            filling = np.trapz(summed_dos_ads[filled_indices],
                               energies_ads[filled_indices]) / np.trapz(
                                   summed_dos_ads, energies_ads)
            print(metal, value, vacancy, filling)
            axp[i].fill_between(summed_dos_ads[filled_indices],
                                energies_ads[filled_indices],
                                alpha=0.25,
                                color='tab:green')
            axp[i].plot(summed_dos_slab,
                        energies_slab,
                        color='tab:purple',
                        alpha=0.3)
            # axp[i].fill_between(summed_dos_slab, energies_slab, color='tab:purple', alpha=0.25)
            axp[i].annotate(r'%sN$_{4} - \left ( \mathregular{d} \right )$' %
                            (metal),
                            xy=(0.3, 0.01),
                            xycoords='axes fraction',
                            fontsize=14,
                            color='tab:purple')
            axp[i].set_ylim([-5, 4])
            # axp[i].set_yticks([])
            axp[i].set_xticks([])
            axp[i].set_xlim([0.0, 1.2])
            if i != 0:
                axp[i].set_yticks([])
            else:
                axp[i].set_ylabel(r'$\epsilon - \epsilon_{f}$ / eV')
                # axp[i].plot([],[], color='tab:blue', lw=10, label=r'CO$_{2}*$')
                axp[i].annotate(r'$\mathregular{CO}_{2}^{*} \left ( \mathregular{s,p} \right)$', xy=(0.2,0.9),\
                     color='tab:green', xycoords='axes fraction', fontsize=14)
                axp[i].legend(loc='best', frameon=False, fontsize=14)
            i += 1

    ## plot the cdd
    cdd = mpimg.imread('schematic/schematic.png')
    axc.imshow(cdd)
    axc.set_xticks([])
    axc.set_yticks([])
    axc.plot([], [], color='b', label=r'$\rho = 0.011 \mathregular{e}$')
    axc.plot([], [], color='orange', label=r'$\rho = -0.011 \mathregular{e}$')
    axc.legend(loc='lower left', frameon=False, fontsize=14)
    # with open('inputs/cdd_Au_211.json', 'r') as handle:
    #     data = json.load(handle)
    # atoms = read('inputs/Au211.traj')
    # axt = axd.twinx()
    # for atom in atoms:
    #     color = jmol_colors[atom.number]
    #     radius = radii[atom.number]
    #     if atom.number in [6,8]:
    #         alpha=1
    #     else:
    #         alpha=0.25
    #     circle = Circle((atom.x, atom.z), radius, facecolor=color,
    #                             edgecolor='k', linewidth=0, alpha=alpha)
    #     # axd.add_patch(circle)
    #     axt.add_patch(circle)
    # axd.plot( data['cdd'], data['z'], '-', color='tab:blue', lw=4)
    # intergral = [np.trapz(data['cdd'][:i], data['z'][:i]) for i in range(len(data['cdd']))]
    # axd.plot(data['z'], intergral, color='tab:red', lw=4)
    # axd.set_xticks([])
    # axd.set_yticks([])
    # axd.set_title(r'$\Delta \rho$')
    # axd.set_xlabel(r'z / $\AA$')
    # axt.set_ylim([-8,10])
    # axt.set_yticks([])

    # pos = axd.imshow(data['xy'])
    # fig.colorbar(pos, ax=axd)
    # axd.set_yticks([])
    # axd.set_xticks([])
    # axd.set_ylabel(r'$\Delta \rho$')
    # axd.set_xlabel(r'z / $\AA$')

    alphabet = list(string.ascii_lowercase)
    ax_all = [ax, axc, axp[0]]  #+ axp
    for i, a in enumerate(ax_all):
        a.annotate(alphabet[i] + ')',
                   xy=(0, 1.1),
                   xycoords='axes fraction',
                   fontsize=20)

    fig.savefig('output/dipole.pdf')
    fig.savefig('output/dipole.png', dpi=300)


if __name__ == '__main__':
    get_plot_params()
    create_output_directory()
    main()
