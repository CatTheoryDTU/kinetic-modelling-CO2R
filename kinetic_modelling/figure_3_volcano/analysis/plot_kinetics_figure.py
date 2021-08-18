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


def get_electronic_energy_points(energy_file, surfaces, facet):
    data = [a.split('\t') for a in energy_file.split('\n')]
    results = {}
    for i, dat in enumerate(data):
        if i == 0:
            continue
        if dat[0] in surfaces and dat[1] == facet:
            dE = float(dat[3])
            results.setdefault(dat[0], {})[dat[2] + '_s'] = dE
        if dat[0] in ['Ni', 'Fe']:
            dE = float(dat[3])
            results.setdefault(dat[0] + '_' + dat[1], {})[dat[2] + '_s'] = dE

    return results

def plot_map(fig, ax, maps, descriptors, points, potential, pH, \
            plot_single_atom, plot_metal, coverage_index=-1, min_val=1e-20, log_scale=True,\
                cmapname='Blues_r', annotate_rate_limiting=False,\
                    coverage_plot=False, plot_cmap=True, inten=1, plot_legend=False):
    """Generate the main kinetic plot; this function keeps getting called for each axis

    :param fig: matplotlib figure
    :type fig: fig
    :param ax: axis to plot on
    :type ax: ax
    :param maps: map to plot on axis
    :type maps: list
    :param descriptors: descriptor names
    :type descriptors: list
    :param points: data points to be plotted on graph
    :type points: dict
    :param potential: SHE potential
    :type potential: float
    :param pH: pH that was used in the run
    :type pH: float
    :param plot_metal: Should the metal points be plotted?
    :type plot_metal: bool
    :param coverage_index: index of the species used to plot the coverage - same indexing as CatMAP, defaults to -1
    :type coverage_index: int, optional
    :param min_val: minimum value to be used, defaults to 1e-20
    :type min_val: float, optional
    :param log_scale: Should a log scale be used?, defaults to True
    :type log_scale: bool, optional
    :param annotate_rate_limiting: should annotations be used?, defaults to False
    :type annotate_rate_limiting: bool, optional
    :param plot_cmap: Plot the colormap?, defaults to True
    :type plot_cmap: bool, optional
    :param inten: Intensity of plot - alpha setting in matplotlib, defaults to 1
    :type inten: int, optional
    """

    ## CHE correction for COOH*
    species = [r'CO_{2}', 'COOH', 'CO']
    CHE_COOH = potential + 0.059 * pH

    x = []
    y = []
    z = []
    z_cov = []  ## for the coverage plot
    for row in maps:
        descriptor, production = row
        dGa, dGb = descriptor
        x.append(dGa)
        y.append(dGb)
        if not coverage_plot:
            rate_CO = np.max(production)
            z.append(rate_CO)
        if coverage_plot:
            z.append(production[coverage_index])

    x = [a + CHE_COOH for a in x]
    z = np.array(z).transpose()
    x = np.array(x)
    y = np.array(y)

    ## Add free energy contributions
    x = x + 1.308
    y = y + 0.766

    # Remove any rates that are very low
    if plot_cmap:
        z1 = [min_val if a_ < min_val else a_ for a_ in z]
        if log_scale:
            x_dense_val = np.linspace(min(x), max(x))
            y_dense_val = np.linspace(min(y), max(y))
            x_dense, y_dense = np.meshgrid(x_dense_val, y_dense_val)
            z_smooth = Rbf(x, y, np.log10(z1), function='linear', smooth=1)
            z_smooth_dense = z_smooth(x_dense, y_dense)
            tcf = ax.contourf(x_dense,
                              y_dense,
                              10**z_smooth_dense,
                              cmap=cmapname,
                              locator=ticker.LogLocator())
        else:
            x_dense_val = np.linspace(min(x), max(x))
            y_dense_val = np.linspace(min(y), max(y))
            x_dense, y_dense = np.meshgrid(x_dense_val, y_dense_val)
            z_smooth = Rbf(x, y, z, function='linear', smooth=1)
            z_smooth_dense = z_smooth(x_dense, y_dense)
            z_smooth_dense = z_smooth_dense.clip(min=0, max=1.)
            tcf = ax.contourf(x_dense,
                              y_dense,
                              z_smooth_dense,
                              cmap=cmapname,
                              levels=np.arange(0, 1.1, 0.1))

        cb = fig.colorbar(tcf, ax=ax)
        if coverage_plot:
            cb.ax.set_title(r'$\theta_{\mathregular{%s}}$ / ML ' %
                            species[coverage_index])
        else:
            cb.ax.set_title(r'TOF / $\mathregular{s}^{-1}$')
        ax.plot([], [], 'o', color='r', label='Fe(SV)')
        ax.plot([], [], 'o', color='g', label='Fe(DV)')
        ax.plot([], [], 'o', color='b', label='Ni(SV)')
        ax.plot([], [], 'o', color='tab:brown', label='Ni(DV)')

        if plot_legend:
            ax.legend(loc='lower right', frameon=False, fontsize=16)
        ax.set_xlabel(r'$\Delta G_{\mathregular{' +
                      descriptors[0].replace('_s', '').replace('2', '_{2}') +
                      '}}(%1.1f\mathregular{V}_{\mathregular{SHE}})$ / eV' %
                      potential)
        ax.set_ylabel(r'$\Delta G_{\mathregular{' +
                      descriptors[1].replace('_s', '').replace('2', '_{2}') +
                      '}}(%1.1f\mathregular{V}_{\mathregular{SHE}})$ / eV' %
                      potential)

    xbounds = ax.get_xbound()
    ybounds = ax.get_ybound()

    all_Ga = []
    all_Gb = []
    for metal in points:
        try:
            color = jmol_colors[atomic_numbers[metal]]
            color = 'k'
        except KeyError:
            color = jmol_colors[atomic_numbers[metal.split('_')[0]]]
        try:
            if points[metal][descriptors[0]] > 1.5:
                continue
            # Elements to scale with
            x = points[metal][descriptors[0]] + CHE_COOH
            y = points[metal][descriptors[1]]

            x = x + 1.308
            y = y + 0.772

            if 'Ni' in metal and plot_single_atom == False:
                continue
            if 'Fe' in metal and plot_single_atom == False:
                continue
            if 'Ni' in metal or 'Fe' in metal and plot_single_atom == True:
                label_info = metal.split('_')
                if int(label_info[1]) == 1 and 'Ni' in metal: color_vac = 'b'
                elif int(label_info[1]) == 2 and 'Ni' in metal:
                    color_vac = 'tab:brown'
                elif int(label_info[1]) == 1 and 'Fe' in metal:
                    color_vac = 'r'
                elif int(label_info[1]) == 2 and 'Fe' in metal:
                    color_vac = 'g'
                ax.plot(
                    x,
                    y,
                    marker='o',
                    markersize=12,
                    color=color_vac,
                )
            elif plot_metal:
                ax.plot(x, y, 'o', markersize=12, color=color)

        except KeyError:
            continue

        if plot_metal:
            if 'Ni' in metal or 'Fe' in metal:
                # continue
                met, vac, dop = metal.split('_')
                ax.annotate('$\mathregular{%sN_{%s}}$' % (met, dop),
                            fontsize=12,
                            color='k',
                            xy=(x + 0.05, y + 0.05),
                            alpha=0.5).draggable()
            else:
                met = metal
                ax.annotate(metal,
                            xy=(x + 0.05, y + 0.05),
                            color=color,
                            fontsize=14).draggable()

        if 'Ni' not in metal and 'Fe' not in metal:
            all_Ga.append(x)
            all_Gb.append(y)

    fit = get_fit_from_points(all_Ga, all_Gb, 1)
    print(fit)
    bounds = ax.get_xbound()

    if plot_metal:
        ax.plot(bounds, fit['p'](bounds), color='k', ls='--', lw=3)
        ax.fill_between(bounds,
                        fit['p'](bounds) + 0.1,
                        fit['p'](bounds) - 0.1,
                        color='k',
                        alpha=0.1)

    ax.plot(bounds, bounds, color='k', ls='-')
    ax.fill_between(bounds,
                    np.array(bounds) + 0.1,
                    np.array(bounds) - 0.1,
                    color='k',
                    alpha=0.1)

    if annotate_rate_limiting:
        ax.arrow(-1.2, 0.5, 0, 1, width=0.05, color='white')
        ax.arrow(0, -1., 1, 0, width=0.05, color='white')
        ax.annotate(
            r'$\mathregular{CO_2}^* \to \mathregular{COOH}^* $ limited',
            xy=(0.4, 0.05),
            xycoords='axes fraction',
            color='white')
        ax.annotate(
            r'$ \mathregular{CO_2}_{(\mathregular{g})} \to \mathregular{CO_2}^*$ limited',
            xy=(0.1, 0.9),
            xycoords='axes fraction',
            color='white')
        ax.annotate('Parity Line',
                    xy=(0.4, 0.2),
                    xycoords='axes fraction',
                    rotation=37,
                    color='k')
        ax.annotate('(211) Scaling',
                    xy=(0.1, 0.4),
                    xycoords='axes fraction',
                    rotation=37,
                    color='k')
    ax.set_xlim(xbounds)
    ax.set_ylim(ybounds)
    # ax.set_xticks(np.arange(-1.5,2,0.5))
    # ax.set_yticks(np.arange(-1.5,2,0.5))


@click.command()
@click.option('--kfiles',
              type=str,
              default='aiida_output/kinetic_model_data.json')
@click.option('--kineticspk', type=str, default='12384')
@click.option('--coveragepk', type=str, default='12251')
def main(kfiles, kineticspk, coveragepk):

    with open(kfiles, 'r') as handle:
        data_tot = json.load(handle)
    data = data_tot[kineticspk]
    # dataco2 = data_tot[coveragepk]

    fig = plt.figure(constrained_layout=True, figsize=(8, 9))
    figs = plt.figure(constrained_layout=True, figsize=(10, 9))
    figp = plt.figure(constrained_layout=True, figsize=(10, 4))
    figt = plt.figure(constrained_layout=True, figsize=(8, 6))
    gs = fig.add_gridspec(6, 1)
    gsS = figs.add_gridspec(2, 2)
    gsP = figp.add_gridspec(1, 2)
    gsT = figt.add_gridspec(1, 1)
    axc = [fig.add_subplot(gs[0:3, 0]),
           fig.add_subplot(gs[3:, 0])]  # computational figure
    axt = figt.add_subplot(gsT[0, 0])  # TPD figure
    axp = [figp.add_subplot(gsP[0, 0]),
           figp.add_subplot(gsP[0, 1])]  # CO poisoning figure
    axd = [figs.add_subplot(gsS[0, 0]),
           figs.add_subplot(gsS[0, 1])]  # DEMS plot with currents
    axco2 = figs.add_subplot(gsS[1, :])  # co2 poisoning plot
    ax = axc  #+  axp #+ axd

    ## open excel books
    dems = xlrd.open_workbook('experiments/DEMS.xlsx')
    tpd = xlrd.open_workbook('experiments/TPD.xlsx')
    cv = xlrd.open_workbook('experiments/CV.xlsx')

    ## plot the TPD graph
    colors = ['tab:red', 'tab:blue', 'tab:green']
    sheet_names_tpd = tpd.sheet_names()
    for i, label in enumerate(sheet_names_tpd):
        sheet = tpd.sheet_by_index(i)
        temp_sheet = sheet.col_slice(0, start_rowx=2, end_rowx=None)
        signal_sheet = sheet.col_slice(-1, start_rowx=2, end_rowx=None)
        temperature = []
        signal = []
        for j in range(len(temp_sheet)):
            temperature.append(float(temp_sheet[j].value) + 273.15)
            signal.append(float(signal_sheet[j].value))
        index = [a for a in range(len(temperature)) if temperature[a] > 300]
        temperature = np.array(temperature)
        signal = np.array(signal)
        axt.plot(temperature[index],
                 signal[index],
                 '-',
                 lw=4,
                 color=colors[i],
                 label=label)
    axt.set_ylabel(r'Signal / arb. units')
    axt.set_yticks([])
    axt.set_xlabel(r'Temperature / K')
    axt.legend(loc='best', frameon=False)
    # axt.annotate('First Peak', color='tab:red', xy=(0.15,0.8), xycoords='axes fraction')
    # axt.annotate('First Peak', color='tab:blue', xy=(0.2,0.5), xycoords='axes fraction')
    # axt.annotate('Second Peak', color='tab:red', xy=(0.4,0.2), xycoords='axes fraction')

    ## plot the CVs
    sheet_names_cv = cv.sheet_names()
    for i, label in enumerate(sheet_names_cv):
        sheet = cv.sheet_by_index(i)
        E_1 = sheet.col_slice(0, start_rowx=2, end_rowx=None)
        cv_1 = sheet.col_slice(1, start_rowx=2, end_rowx=None)
        E_2 = sheet.col_slice(2, start_rowx=2, end_rowx=None)
        cv_2 = sheet.col_slice(3, start_rowx=2, end_rowx=None)
        E_3 = sheet.col_slice(4, start_rowx=2, end_rowx=None)
        cv_3 = sheet.col_slice(5, start_rowx=2, end_rowx=None)
        datap = []
        for j in range(len(E_1)):
            datap.append([
                E_1[j].value, cv_1[j].value, E_2[j].value, \
                cv_2[j].value, E_3[j].value, cv_3[j].value
            ])
        datap = np.array(datap)
        e_N2, cv_N2, e_H2, cv_H2, e_CO, cv_CO = datap.transpose()
        axp[i].plot(e_N2, cv_N2, '-', lw=4, color='tab:blue', label=r'N$_{2}$')
        axp[i].plot(e_H2, cv_H2, '-', lw=4, color='tab:red', label=r'H$_{2}$')
        axp[i].plot(e_CO, cv_CO, '-', lw=4, color='tab:green', label=r'CO')

        axp[i].set_ylabel(r'j$_{\mathregular{tot}}$ / mAcm$^{-2}$')
        axp[i].set_xlabel(r'Potential / V vs. RHE')
        if i == 0:
            axp[i].legend(loc='best', frameon=False)
        axp[i].set_title(label)

    ## DEMS plot
    sheet_names_dems = dems.sheet_names()
    for i, label in enumerate(sheet_names_dems):
        sheet = dems.sheet_by_index(i)
        E_sheet = sheet.col_slice(0, start_rowx=2, end_rowx=None)
        j_sheet = sheet.col_slice(1, start_rowx=2, end_rowx=None)
        mz_sheet = sheet.col_slice(2, start_rowx=2, end_rowx=None)
        datap = []
        for j in range(len(E_sheet)):
            datap.append([
                E_sheet[j].value,
                j_sheet[j].value,
                mz_sheet[j].value,
            ])
        datap = np.array(datap)
        E, current, mz = datap.transpose()
        axd[i].plot(E, mz, color='tab:red')
        ax2d = axd[i].twinx()
        ax2d.plot(E, current, lw=4, color='tab:green')
        axd[i].set_ylabel(r'Normalised Intensity / arb. unit.',
                          color='tab:red')
        ax2d.set_ylabel(r'j$_{\mathregular{geo}}$ / mAcm$^{-2}$',
                        color='tab:green')
        axd[i].set_xlabel(r'Potential / V vs. RHE')
        axd[i].set_title(label)

        axd[i].arrow(-0.4,
                     0,
                     0.4,
                     0.,
                     head_width=0.01,
                     head_length=0.05,
                     color='tab:red')

    ## computational plot
    data_points = get_electronic_energy_points(data['energy_file'],
                                               data['surfaces'],
                                               data['facet'][0])
    # data_points_co2 = get_electronic_energy_points(dataco2['energy_file'], dataco2['surfaces'], dataco2['facet'][0])
    ## set limits for axc
    for a in axc:
        a.set_xlim([-1.5, 2])
        a.set_ylim([-1.5, 2])

    plot_map(
        fig=fig,
        ax=axc[0],
        maps=data['production_rate'],
        descriptors=data['descriptors'],
        points=data_points,
        potential=data['potential'],
        pH=data['pH'],
        plot_single_atom=True,
        cmapname='coolwarm',
        plot_metal=True,
        annotate_rate_limiting=True,
    )
    plot_map(
        fig=fig,
        ax=axc[1],
        maps=data['coverage_map'],
        descriptors=data['descriptors'],
        points=data_points,
        potential=data['potential'],
        pH=data['pH'],
        plot_single_atom=True,
        cmapname='Oranges',
        plot_metal=True,
        log_scale=False,
        coverage_plot=True,
        coverage_index=-1,
        plot_legend=True,
    )
    axc[1].annotate(r'CO$^*$ poisoned',
                    xy=(0.03, 0.9),
                    color='white',
                    xycoords='axes fraction')

    # plot_map(
    #     fig=fig,
    #     ax=axco2,
    #     maps=dataco2['coverage_map'],
    #     descriptors=dataco2['descriptors'],
    #     points=data_points_co2,
    #     potential=dataco2['potential'],
    #     pH=dataco2['pH'],
    #     plot_single_atom=True,
    #     cmapname='Oranges',
    #     plot_metal=True,
    #     log_scale=False,
    #     coverage_plot=True,
    #     coverage_index=0,
    #     plot_legend=True,
    # )
    # axco2.annotate(r'CO$_{2}^*$ poisoned', xy=(0.6, 0.2), color='white', xycoords='axes fraction')
    # axco2.set_ylim([-1,1])
    # axco2.set_xlim([-1,1])

    ## Label the diagram
    alphabet = list(string.ascii_lowercase)
    for i, a in enumerate(ax):
        a.annotate(alphabet[i] + ')',
                   xy=(0.0, 1.1),
                   xycoords='axes fraction',
                   fontsize=20)
    for i, a in enumerate(axd + [axco2]):
        a.annotate(alphabet[i] + ')',
                   xy=(0.0, 1.1),
                   xycoords='axes fraction',
                   fontsize=20)

    fig.savefig('output_figure/figure_kinetics.pdf')
    figs.savefig('output_figure/SI_figure_CO2.pdf')
    figp.savefig('output_figure/SI_figure_CO_poison.pdf')
    figt.savefig('output_figure/TPD.pdf')
    plt.close(figs)
    plt.close(figp)
    plt.close(figt)
    plt.show()


if __name__ == '__main__':
    create_output_directory('output_figure')
    get_plot_params()
    main()
