# -*- coding: utf-8 -*-
# pylint: disable-all
import os.path as op
import matplotlib.pyplot as plt
from ase.calculators.vasp import VaspDos
from plot_params import get_plot_params
get_plot_params()

if __name__ == '__main__':

    folders = {
        1: 'inputs/vacancy_1',
        2: 'inputs/vacancy_2',
        #        3:'inputs/vacancy_3',
        4: 'inputs/vacancy_4',
    }
    fig, ax = plt.subplots(1, 1, figsize=(4, 6), constrained_layout=True)

    for system, folder in folders.items():
        dos = VaspDos(op.join(folder, 'DOSCAR'))
        # get Fermi level from DOSCAR
        dosfile = open(op.join(folder, 'DOSCAR'))
        for i, line in enumerate(dosfile):
            if i == 5:
                values = line.split()
                break

        # get energies and subtract Fermi level
        efermi = float(values[-2])

        energy = dos.energy - efermi
        # get all DOS
        total_dos = dos.dos

        ax.plot(energy, total_dos, lw=3)

    ax.set_xlim([-2, 2])
    ax.set_ylim([-0.2, 5])
    #     ax.set_yticks([])
    ax.set_ylabel(r'Density of States / a.u.')
    ax.set_xlabel(r'Energy / eV')
    fig.savefig('output/dos.png')
