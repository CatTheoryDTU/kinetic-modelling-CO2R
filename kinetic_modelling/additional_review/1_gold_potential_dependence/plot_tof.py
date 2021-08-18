# -*- coding: utf-8 -*-
# pylint: disable-all
# pylint: skip-file
import json
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np
from plot_params import get_plot_params
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
get_plot_params()


def get_tafel_slope(potential, tof, range_val=[]):
    potential = np.array(potential)
    tof = np.array(tof)
    if range_val:
        index_potential = [
            a for a in range(len(potential))
            if np.min(range_val) < potential[a] < np.max(range_val)
        ]
        potential = potential[index_potential]
        tof = tof[index_potential]

    fit = np.polyfit(np.log10(tof), potential, 1)
    return -1 * fit[0] * 1e3


if __name__ == '__main__':

    data = json.load(
        open(f"output/node_{15157}_surface_{'Au'}_facet_{'211'}.json"))
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)

    ## plot the computation tof
    production_rate = data
    potential_ph = [production_rate[i][0] for i in range(len(production_rate))]
    rate = [production_rate[i][1][1] for i in range(len(production_rate))]
    tof_comp = np.array(rate)
    she_potential, pH = np.array(potential_ph).T

    ax.plot(she_potential, tof_comp, '-')
    tafel_slope = get_tafel_slope(she_potential,
                                  tof_comp,
                                  range_val=[-0.6, -1.])
    ax.set_ylabel(r'TOF / s$^{-1}$')
    ax.set_xlabel(r'Potential / V vs. SHE')

    ax.set_yscale('log')
    ax.yaxis.grid(True, which='minor')
    ax.set_title('Au(211)')
    ax.annotate('%d mV/dec' % tafel_slope,
                xy=(0.4, 0.7),
                xycoords='axes fraction')

    fig.savefig('output/tof_computations.png')
