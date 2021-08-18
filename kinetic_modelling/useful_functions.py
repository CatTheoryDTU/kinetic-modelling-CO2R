# -*- coding: utf-8 -*-
# pylint: disable-all
def get_fit_from_points(x, y, order):
    import numpy as np
    fit = np.polyfit(x, y, order)
    p = np.poly1d(fit)
    return {'fit': fit, 'p': p}


def get_vasp_nelect0(atoms):
    import pickle
    import numpy as np
    with open('../utilities/nelect0.pickle', 'rb') as handle:
        nelect0 = pickle.load(handle)
    default_nelect = []
    for i in range(len(atoms)):
        default_nelect.append(nelect0[atoms[i].symbol])
    default_nelect = np.array(default_nelect)
    return default_nelect.sum()


def get_reference_energies(database_file):
    results = {}
    from ase.db import connect
    database = connect(database_file)
    for row in database.select():
        state = row.states
        pw = row.pw
        functional = row.functional
        results.setdefault(state, {}).setdefault(functional, {}).setdefault(
            pw, {})['energy'] = row.energy
        results[state][functional][pw]['vibrations'] = row.data.vibrations
        results[state][functional][pw]['atoms'] = row.toatoms()

    return results
