# -*- coding: utf-8 -*-
# pylint: disable-all
import os.path as op
import matplotlib.pyplot as plt
from ase.calculators.vasp import VaspDos
from plot_params import get_plot_params
import json
import numpy as np
get_plot_params()

if __name__ == '__main__':

    folders = {
        'Pt': 'cdd_Pt_211.json',
        'Au': 'cdd_Au_211.json',
        'FeN4': 'Fe_cdd.json',
    }
    system_name = {'Pt': 'Pt(211)', 'Au': 'Au(211)', 'FeN4': r'FeN$_{4}$'}
    fig, ax = plt.subplots(len(folders) + 1,
                           1,
                           figsize=(7, 10),
                           constrained_layout=True)

    for i, (system, folder) in enumerate(folders.items()):
        data = json.load(open(op.join('inputs', folder)))
        xy_av_data = data['xy']
        z_av_cdd = data['cdd']
        z_coord = data['z']
        # ax[i].imshow(xy_av_data, origin='lower', extent=(0,20,0,10))
        if i == 0:
            axt = ax[i].imshow(xy_av_data)
            fig.colorbar(axt, ax=ax[i])
            ax[i].set_xticks([])
            ax[i].set_yticks([])
            ax[i].annotate(r'$\mathregular{CO}_{2}^*$ on Pt(211)',
                           xy=(0.1, 0.8),
                           xycoords='axes fraction',
                           color='white')
            ax[i].set_ylabel(r'$\rho_{z}$ / $e/\AA^3$')
        ax[i + 1].plot(
            z_coord,
            z_av_cdd,
        )
        integrated = np.trapz(z_av_cdd, z_coord)
        print(integrated)
        ax[i + 1].annotate(system_name[system],
                           xy=(0.1, 0.8),
                           xycoords='axes fraction',
                           color='tab:blue')
        # ax[i+1].legend(loc='best')
        ax[i + 1].set_ylabel(r'$\rho_{xy}(z)$ /  $e/\AA^3$')
        ax[i + 1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax[-1].set_xlabel(r'$L_{z}$ / $\AA$')
    fig.savefig('output/density.png')
