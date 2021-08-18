# -*- coding: utf-8 -*-
# pylint: disable-all
import argparse
from pathlib import Path
from ase.io import read
from ase.db import connect
from pprint import pprint
import json
from useful_classes import bcolors
import os


def main(database, fields_to_choose, displacement, direction, db_params):

    results = {}

    # atomsFS = read(os.path.join('sampling_%s'%metal, 'facet_%s'%facet, 'state_implicit_CO2', 'charge_2.00', 'CONTCAR'))
    # atomsIS = read(os.path.join('sampling_%s'%metal, 'facet_%s'%facet, 'state_CO2_gas', 'CONTCAR'))

    method = ForceExtrapolation(
        findiffdb=database,
        fields_to_choose=fields_to_choose,
        displacement=displacement,
        direction=direction,
        db_params=db_params,
        atomsIS=atomsIS,
        atomsFS=atomsFS,
    )

    method.get_dFdG()
    print(f'{bcolors.OKBLUE} dFdG {metal} {bcolors.ENDC}')
    pprint(method.dFdG)

    method.get_dmudR()
    # print(f'{bcolors.OKBLUE} dmudR {bcolors.ENDC}')
    # pprint(method.dmudR)

    method.get_reaction_path()
    method.get_q()

    print(f'{bcolors.OKBLUE} modes {metal} {bcolors.ENDC}')
    pprint(method.modes)
    print(f'{bcolors.OKBLUE} q {metal} {bcolors.ENDC}')
    pprint(method.q)

    results.setdefault(facet, {}).setdefault(
        metal, {})['q_eff'] = method.q['state_implicit_CO2'][0]

    with open(os.path.join(outdir, 'q_exp.json'), 'w') as handle:
        json.dump(results, handle)
