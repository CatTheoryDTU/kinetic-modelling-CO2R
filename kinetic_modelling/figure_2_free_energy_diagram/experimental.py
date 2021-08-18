# -*- coding: utf-8 -*-
# pylint: disable-all
"""Plots experimental data for the rates against potential and pH.
"""
import xlrd
import matplotlib
import numpy as np
from kinetic_modelling.useful_functions import get_fit_from_points  # pylint: disable=no-name-in-module, import-error


def plot_experimental_data(  # pylint: disable=too-many-arguments, too-many-locals, too-many-statements, too-many-branches, invalid-name
        f,
        ax,
        annotation,
        title,
        color,
        fit_min=-0.75,
        fit_lim=-1.,
        fit_all=True,
        plot_line=True,
        marker='o',
        pH_material=None):
    """Plots experimental data for the rates against potential and pH."""

    wb = xlrd.open_workbook(f)  # pylint: disable=invalid-name

    sheet_names = wb.sheet_names()
    cmap = matplotlib.cm.get_cmap('coolwarm', len(sheet_names))

    fit_potential_all = []
    fit_log_current_all = []
    fit_current_all = []

    tafel_all = []
    for i, label in enumerate(sheet_names):
        sheet = wb.sheet_by_index(i)
        pot_sheet = sheet.col_slice(0, start_rowx=2, end_rowx=None)
        current_sheet = sheet.col_slice(-1, start_rowx=2, end_rowx=None)

        potential = []
        current = []

        for j in range(len(pot_sheet)):  # pylint: disable=consider-using-enumerate
            try:
                potential.append(float(pot_sheet[j].value))
            except ValueError:
                continue
            current.append(float(current_sheet[j].value))


        ax.plot(potential, current, marker, \
                markersize=5,
                label=label.replace('KHCO3', '').replace('Phosphate', ''), color=cmap(i))
        potential = np.array(potential)
        current = np.array(current)

        # Fit the Tafel slope
        fit_index = [
            i for i in range(len(potential))
            if fit_min > potential[i] > fit_lim
        ]
        fit_potential = potential[fit_index]
        fit_log_current = np.log10(current[fit_index])
        fit_current = current[fit_index]

        if fit_all:
            fit_potential_all += fit_potential.tolist()
            fit_log_current_all += fit_log_current.tolist()
            fit_current_all += fit_current.tolist()

        if not fit_all:
            if 'pH 6.4' in label.replace('KHCO3', '').replace('Phosphate', ''):
                if pH_material == 'Co':
                    fit_index = [
                        i for i in range(len(potential))
                        if -0.5 > potential[i] > -0.9
                    ]
                elif pH_material == 'Ni':
                    fit_index = [
                        i for i in range(len(potential))
                        if -0.88 > potential[i] > -1.1
                    ]

                fit_potential = potential[fit_index]
                fit_log_current = np.log10(current[fit_index])
                fit_current = current[fit_index]
            if 'pH 6.8' in label.replace('KHCO3', '').replace('Phosphate', ''):
                if pH_material == 'Co':
                    fit_index = [
                        i for i in range(len(potential))
                        if -0.5 > potential[i] > -1.
                    ]
                elif pH_material == 'Ni':
                    fit_index = [
                        i for i in range(len(potential))
                        if -0.85 > potential[i] > -1.1
                    ]
                fit_potential = potential[fit_index]
                fit_log_current = np.log10(current[fit_index])
                fit_current = current[fit_index]
            fit_tafel = get_fit_from_points(fit_log_current, fit_potential, 1)
            fit_plot = get_fit_from_points(fit_potential, fit_log_current, 1)
            ax.plot(potential,
                    10**fit_plot['p'](potential),
                    color=cmap(i),
                    alpha=0.5)
            tafel_all.append(-1 * fit_tafel['fit'][0] * 1000)

    if not fit_all:
        min_tafel = min(tafel_all)
        max_tafel = max(tafel_all)
        ax.annotate(r'$ %1.0f - %1.0f \frac{mV}{dec}$'%(min_tafel, max_tafel), \
            xy=(0.55, 0.85), xycoords='axes fraction', color='k',)
    ax.set_yscale('log')
    ax.set_ylabel(r'$j_{\mathregular{CO}}$ / mAcm$^{-2}$')
    ax.set_xlabel(r'Potential / V vs NHE')
    if annotation == 'pH dependent':
        ax.annotate(annotation,
                    xy=(0.8, 0.4),
                    xytext=(0.05, 0.4),
                    xycoords='axes fraction',
                    arrowprops={
                        'arrowstyle': '<->',
                        'color': 'k',
                        'lw': 2
                    },
                    va='center',
                    fontsize=14)
    else:
        if plot_line:
            ax.annotate(annotation,
                        xy=(0.05, 0.4),
                        xycoords='axes fraction',
                        fontsize=14)
    ax.set_xlim([-1.5, -0.4])
    ax.set_ylim([1e-3, 1e3])

    ax.tick_params(direction='out',
                   length=6,
                   width=2,
                   colors='k',
                   grid_color='r',
                   grid_alpha=0.5)

    ax.legend(loc='lower left', frameon=False, fontsize=14)

    if fit_all:
        fit_plot = get_fit_from_points(fit_potential_all, fit_log_current_all,
                                       1)
        fit_tafel = get_fit_from_points(fit_log_current_all, fit_potential_all,
                                        1)
        plot_potential = np.linspace(-0.5, -1.3)
        if plot_line:
            ax.plot(plot_potential,
                    10**fit_plot['p'](plot_potential),
                    ls='-',
                    color='k',
                    lw=4)
            tafel_slope = -1 * fit_tafel['fit'][0] * 1000
            ax.annotate(r'$%1.0f \frac{\mathregular{mV}}{\mathregular{dec}}$' %
                        tafel_slope,
                        xy=(0.2, 0.82),
                        xycoords='axes fraction',
                        color='k')

    ax.annotate(title, color=color, xy=(0.4, 0.9), xycoords='axes fraction')
