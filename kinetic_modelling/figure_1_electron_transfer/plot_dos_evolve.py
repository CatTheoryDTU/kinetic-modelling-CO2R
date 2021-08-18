# -*- coding: utf-8 -*-
# pylint: skip-file
"""
Generate Figure 1 of the paper.
"""
import string  # pylint: disable=syntax-error
from pathlib import Path
import numpy as np
from ase.db import connect
from ase.data.colors import jmol_colors
from ase.data import covalent_radii as radii
import click
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from kinetic_modelling.plot_params import get_plot_params_Arial as get_plot_params  # pylint: disable=import-error, no-name-in-module


def parsedb(database, results):
    """Parse the database to get the density of states data."""
    for row in database.select():
        results.setdefault(row.states.replace('state_', ''),
                           {})['pdos'] = row.data.pdos.pdos
        results.setdefault(row.states.replace('state_', ''),
                           {})['energy'] = row.data.pdos.energies
        results.setdefault(row.states.replace('state_', ''),
                           {})['integrated_dos'] = row.data.pdos.integrated_dos
        results.setdefault(row.states.replace('state_', ''),
                           {})['total_dos'] = row.data.pdos.total_dos
        results.setdefault(row.states.replace('state_', ''), {})['wf'] = row.wf
        results.setdefault(row.states.replace('state_', ''),
                           {})['atoms'] = row.toatoms()
        results.setdefault(row.states.replace('state_', ''),
                           {})['magmom'] = row.magmom


def lorentz_dos(a_numer, b1_denom, energy):
    """Create a lorentzian dos."""
    delta = 1 / np.pi * (a_numer / ((energy - b1_denom)**2 + a_numer**2))
    return delta


@click.command()
@click.option('--mncdb',
              default='databases/single_atom_rls.db',
              help='Database for single atoms')
def main(mncdb):  # pylint: disable=too-many-locals, too-many-statements, too-many-branches
    """Main function for generating the plot.

    :param mncdb: Database for the MNC.
    :type mncdb: str
    """

    fig = plt.figure(constrained_layout=True, figsize=(14, 8))
    gs = fig.add_gridspec(2, 10, wspace=0.05)  # pylint: disable=invalid-name
    ax = []  # pylint: disable=invalid-name
    for i in range(1):
        temp_ax = []
        for j in range(5):
            temp_ax.append(fig.add_subplot(gs[i + 1, 2 * j:2 * j + 2]))
        ax.append(temp_ax)
    ax = np.array(ax)  # pylint: disable=invalid-name
    axl = fig.add_subplot(gs[0, 5:])  # pylint: disable=invalid-name
    axr = fig.add_subplot(gs[0, 0:5])  # pylint: disable=invalid-name

    results = {}
    results['mnc'] = {}
    parsedb(connect(mncdb), results['mnc'])

    def plot_lorentz_dos():
        """ Plot the idealised density of states and the rate of electron transfer. """
        energy_range = np.linspace(-3, 3, 500)
        spec_lorentz = [[0.01, 0], [0.1, 0], [0.5, 0], [1, 0], [2, 0]]
        color = [
            'tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple'
        ]
        for m, spec in enumerate(spec_lorentz):  # pylint: disable=invalid-name
            peak = lorentz_dos(*spec, energy_range)
            rate = 2 * np.pi / HBAR * spec[0]
            axl.plot(peak, energy_range, '-', color=color[m])
            axl.fill_between(peak, energy_range, color=color[m], alpha=0.25)
            axr.plot(
                spec[0],
                rate,
                'o',
                markersize=14,
                color=color[m],
            )
        axl.set_ylabel(r'Energy / eV')
        axl.set_xlim([0, 1])
        axr.axhline(1e12, color='k', ls='--')
        axr.annotate('Diffusion rate',
                     xy=(0.5, 0.12),
                     xycoords='axes fraction')
        axr.set_yscale('log')
        axl.set_xticks([])
        axr.set_ylabel(r'Rate / s$^{-1}$')
        axr.set_xlabel(r'Width / eV')
        annotation = r'$\Delta = \Sigma_{k} V_{ak}^2 \delta \left ( \epsilon - \epsilon_{k} \right )$'
        axl.annotate(annotation, xy=(0.4, 0.8), xycoords='axes fraction')

    # Plot the idealised density of states and the rate of electron transfer.
    plot_lorentz_dos()

    plot_states = {}
    plot_states['mnc'] = ['00', '04', '07', '08', '09']
    width = [[-0.1, 0.1], [-0.45, 0.3], [-0.8, 0.3]]

    for i, species in enumerate(results):  # pylint: disable=too-many-nested-blocks
        j = 0
        index = 0  # pylint: disable=invalid-name
        all_f = []
        for state in results[species]:
            energy = np.array(results[species][state]['energy'])
            filled_indices = [
                a for a in range(len(energy)) if energy[a] <= 0.0
            ]
            pdos_co2 = np.array(results[species][state]['pdos']['co2']['+']) \
                    +np.array(results[species][state]['pdos']['co2']['-'])
            co2 = pdos_co2[0].sum(axis=0)
            filling = np.trapz(co2[filled_indices],
                               energy[filled_indices]) / np.trapz(co2, energy)
            all_f.append(filling)
            if state in plot_states[species]:
                ax[i, j].plot(co2, energy, color='tab:blue', alpha=0.5)
                ax[i, j].fill_between(co2, energy, color='tab:blue', alpha=0.1)
                ax[i, j].set_ylim([-8, 5])
                ax[i, j].set_xlim([0.0, 0.1])
                axins = ax[i, j].inset_axes([0.5, 0.5, 0.47, 0.47])
                axins.plot(co2, energy, color='tab:blue')
                axins.fill_between(co2, energy, color='tab:blue', alpha=0.25)
                axins.set_xlim([0.0, 0.006])
                axins.set_ylim([-2, 2])
                mark_inset(ax[i, j],
                           axins,
                           loc1=2,
                           loc2=3,
                           fc='none',
                           lw=3,
                           ec='k')

                axins.set_xticklabels('')
                axins.set_yticklabels('')
                if j == 0:
                    ax[i, j].set_ylabel(r'$\epsilon - \epsilon_{f}$ / eV')
                else:
                    ax[i, j].set_yticks([])
                ax[i, j].set_xticks([])
                ax[i, j].axhline(0, ls='--')
                if j > 1 and species == 'mnc':
                    ax[i, j].axhspan(*width[index], color='tab:red', alpha=0.5)
                    axins.axhspan(*width[index], color='tab:red', alpha=0.5)
                    index += 1

                ## plot the atoms as insets
                atoms = results[species][state]['atoms']
                axins4 = inset_axes(ax[i, j],
                                    width='40%',
                                    height='40%',
                                    loc=4,
                                    borderpad=1)
                for atom in atoms:
                    if species == 'mnc':
                        if atom.index not in [
                                26, 32, 33, 27, 31, 30, 28, 29, 18, 24, 23, 11,
                                17
                        ]:
                            continue

                    color = jmol_colors[atom.number]
                    radius = radii[atom.number]
                    if species == 'mnc':
                        circle = Circle((atom.y, atom.z),
                                        radius,
                                        facecolor=color,
                                        edgecolor='k',
                                        linewidth=1)
                    else:
                        circle = Circle((atom.x, atom.y),
                                        radius,
                                        facecolor=color,
                                        edgecolor='k',
                                        linewidth=1)
                    axins4.add_patch(circle)

                axins4.axis('equal')
                axins4.set_xticks([])
                axins4.set_yticks([])
                axins4.axis('off')

                j += 1

            ax[0, 2].annotate('TS', xy=(0.25, 0.1), xycoords='axes fraction')

    fig.suptitle(
        r'Rate of e$^-$ transfer on MNC $ \approx 10^{14} \mathregular{s}^{-1}$'
    )

    fig.text(0.6, 0.46, '(s,p) PDOS / arb. units.', ha='center')
    alphabet = list(string.ascii_lowercase)
    for i, ax_ in enumerate([axr, axl] + ax.flatten().tolist()):
        ax_.annotate(alphabet[i] + ')',
                     xy=(0.1, 0.9),
                     xycoords='axes fraction',
                     fontsize=18)
    fig.savefig('output/dos.pdf')
    fig.savefig('output/dos.png', dpi=300)


if __name__ == '__main__':
    Path('output').mkdir(exist_ok=True)
    get_plot_params()
    HBAR = 6.582 * 1e-16  # eV.s
    main()  # pylint: disable=no-value-for-parameter
