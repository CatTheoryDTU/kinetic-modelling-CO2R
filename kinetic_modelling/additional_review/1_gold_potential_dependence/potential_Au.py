# -*- coding: utf-8 -*-
# pylint: disable-all
# pylint: skip-file
from os import path
from aiida_catmap import helpers
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
import click
from aiida.orm import SinglefileData, List, Dict, Int, Float, Str

INPUT_DIR = path.join(path.dirname(path.realpath(__file__)), 'input_files')


def run_calculation(catmap_code, energy_pk, label):

    if not catmap_code:
        # get code
        computer = helpers.get_computer()
        catmap_code = helpers.get_code(entry_point='catmap', computer=computer)

    # Prepare input parameters

    energies = load_node(energy_pk)
    species_definitions = {}
    species_definitions['CO2_g'] = {'pressure': 0.2}
    species_definitions['CO_g'] = {'pressure': 0.1}
    species_definitions['H2_g'] = {'pressure': 1}
    species_definitions['H2O_g'] = {'pressure': 0.1}

    ## NiNC
    species_definitions['CO2_s'] = {
        'nsites': 1,
        'sigma_params': [0.01916411, 0.39391097]
    }
    species_definitions['COOH_s'] = {
        'nsites': 1,
        'sigma_params': [0.00425954, -0.19811406]
    }
    species_definitions['CO_s'] = {
        'nsites': 1,
        'sigma_params': [-0.00293112, -0.1824965]
    }
    species_definitions['s'] = {'site_names': ['211'], 'total': 1}
    surfaces = ['Au']
    Upzc = 0.4

    ideal_gas_params = {
        'CO2_g': [2, 'linear', 0],
        'CO_g': [1, 'linear', 0],
        'H2_g': [2, 'linear', 0],
        'H2O_g': [2, 'nonlinear', 0],
    }
    scaling_constraint_dict = {
        'COOH_s': ['+', '+', None],
        'CO2_s': ['+', '+', None],
        'CO_s': ['+', '+', None],
    }

    # set up calculation
    sigma_input = ['CH', 25]  # Helmholtz capacitance

    inputs = {
        'code':
        catmap_code,
        'energies':
        energies,
        'scaler':
        Str('ThermodynamicScaler'),
        'rxn_expressions':
        List(list=[
            'CO2_g + *_s <-> CO2_s',
            'CO2_s + + H_g + ele_g <-> COOH_s',
            'COOH_s + H_g + ele_g <-> CO_s + H2O_g',
            'CO_s <-> CO_g + *_s',
            'H2_g <-> H2_g',
        ]),
        'surface_names':
        List(list=surfaces),
        'descriptor_names':
        List(list=['voltage', 'pH']),
        'descriptor_ranges':
        List(list=[[-1.5, -0.0], [2, 2]]),
        'resolution':
        Int(50),
        'temperature':
        Float(300),
        'species_definitions':
        Dict(dict=species_definitions),
        'gas_thermo_mode':
        Str('ideal_gas'),
        'adsorbate_thermo_mode':
        Str('harmonic_adsorbate'),
        'electrochemical_thermo_mode':
        List(list=['simple_electrochemical', 'surface_charge_density']),
        'decimal_precision':
        Int(150),
        'tolerance':
        Float(1e-20),
        'max_rootfinding_iterations':
        Int(100),
        'max_bisections':
        Int(3),
        'ideal_gas_params':
        Dict(dict=ideal_gas_params),
        'sigma_input':
        List(list=sigma_input),
        'Upzc':
        Float(Upzc),
        'scaling_constraint_dict':
        Dict(dict=scaling_constraint_dict),
        #'numerical_solver':Str('coverages'),
        'metadata': {
            'label': label,
            'description': 'Submissions with AiiDA CatMAP',
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    from aiida.engine import submit
    future = engine.submit(CalculationFactory('catmap'), **inputs)
    # grp = Group.get(label='mkm_doped_graphene/TM_test_calcs')
    # grp = Group.get(label='nipc_stability')
    # grp.add_nodes(future)
    # result = engine.run(CalculationFactory('catmap'), **inputs)

    # computed_result = result['log'].get_content()


@click.command()
@click.option('--energy_pk')
@click.option('--label')
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code, energy_pk, label):
    """Run example.

    Example usage: $ ./example_01.py --code diff@localhost

    Alternative (creates diff@localhost-test code): $ ./example_01.py

    Help: $ ./example_01.py --help
    """
    run_calculation(code, energy_pk, label)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
