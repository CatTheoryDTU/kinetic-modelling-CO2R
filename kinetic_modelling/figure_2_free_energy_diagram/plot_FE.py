# -*- coding: utf-8 -*-
# pylint: disable-all
import xlrd
import matplotlib
import numpy as np
from useful_functions import get_fit_from_points
import matplotlib.pyplot as plt
from useful_functions import create_output_directory
from plot_params import get_plot_params
import string

if __name__ == '__main__':

    create_output_directory()
    get_plot_params()

    f = 'inputs/pH_effect_NiNC.xlsx'
    wb = xlrd.open_workbook(f)

    fig, ax = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

    sheet_names = wb.sheet_names()
    cmap = matplotlib.cm.get_cmap('coolwarm', len(sheet_names))

    for i, label in enumerate(sheet_names):
        sheet = wb.sheet_by_index(i)
        pot_sheet = sheet.col_slice(0, start_rowx=2, end_rowx=None)
        tot_current_sheet = sheet.col_slice(2, start_rowx=2, end_rowx=None)
        fe_h2_sheet = sheet.col_slice(3, start_rowx=2, end_rowx=None)
        fe_co_sheet = sheet.col_slice(4, start_rowx=2, end_rowx=None)
        h2_current_sheet = sheet.col_slice(5, start_rowx=2, end_rowx=None)

        data = []

        for j in range(len(pot_sheet)):
            data.append([
                pot_sheet[j].value,
                tot_current_sheet[j].value,
                fe_h2_sheet[j].value,
                fe_co_sheet[j].value,
                h2_current_sheet[j].value,
            ])

        data = np.array(data)
        u_she, j_tot, fe_h2, fe_co, j_h2 = data.transpose()

        ax[0, 0].plot(u_she,
                      j_tot,
                      'o',
                      markersize=5,
                      label=label.replace('KHCO3',
                                          '').replace('Phosphate', ''),
                      color=cmap(i))
        ax[0, 0].set_ylabel(r'j$_{\mathregular{tot}}$ / mAcm$^{-2}$')
        ax[0, 0].set_xlabel(r'Potential / V vs. NHE')
        ax[0, 0].legend(loc='best', frameon=False, fontsize=12)

        ax[0, 1].plot(u_she,
                      j_h2,
                      'o',
                      markersize=5,
                      label=label.replace('KHCO3',
                                          '').replace('Phosphate', ''),
                      color=cmap(i))
        ax[0, 1].set_ylabel(r'j$_{\mathregular{H}_{2}}$ / mAcm$^{-2}$')
        ax[0, 1].set_xlabel(r'Potential / V vs. NHE')

        ax[1, 0].plot(u_she,
                      fe_h2,
                      'o',
                      markersize=5,
                      label=label.replace('KHCO3',
                                          '').replace('Phosphate', ''),
                      color=cmap(i))
        ax[1, 0].set_ylabel(r'FE$_{\mathregular{H}_{2}}$ / %')
        ax[1, 0].set_xlabel(r'Potential / V vs. NHE')

        ax[1, 1].plot(u_she,
                      fe_co,
                      'o',
                      markersize=5,
                      label=label.replace('KHCO3',
                                          '').replace('Phosphate', ''),
                      color=cmap(i))
        ax[1, 1].set_ylabel(r'FE$_{\mathregular{CO}}$ / %')
        ax[1, 1].set_xlabel(r'Potential / V vs. NHE')

    alphabet = list(string.ascii_lowercase)
    for i, a in enumerate(ax.flatten()):
        a.annotate(alphabet[i] + ')',
                   xy=(0.05, 1.05),
                   xycoords='axes fraction',
                   fontsize=20)

    fig.savefig('output/figure_FE.pdf')
