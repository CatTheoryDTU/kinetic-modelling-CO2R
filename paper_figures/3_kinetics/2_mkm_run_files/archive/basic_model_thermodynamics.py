from os import path
from aiida_catmap import helpers
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
import click
from aiida.orm import SinglefileData, List, Dict, Int, Float, Str

INPUT_DIR = path.join(path.dirname(path.realpath(__file__)), 'input_files')


def run_calculation(catmap_code, potential, pH, energy_pk, label):
    """Run AiiDA CatMAP calculation

    :param catmap_code: Code object for CatMAP
    :type catmap_code: Code
    :param potential: SHE potential
    :type potential: float
    :param pH: pH value
    :type pH: float
    :param energy_pk: pk of the energy file
    :type energy_pk: int
    :param label: Label to denote the computation
    :type label: str
    """
    if not catmap_code:
        # get code
        computer = helpers.get_computer()
        catmap_code = helpers.get_code(entry_point='catmap', computer=computer)

    # Prepare input parameters

    energies = load_node(energy_pk)

    species_definitions = {}
    species_definitions['CO2_g'] = {'pressure':0.2} 
    species_definitions['CO_g'] = {'pressure':0.1}
    species_definitions['H2_g'] = {'pressure':0.1}
    species_definitions['H2O_g'] = {'pressure':0.1}
    species_definitions['s'] = {'site_names': ['211'], 'total':1}

    scaling_constraint_dict = {
                           'CO_s':[None,None,None],
                           'COOH_s':[None,None,None],
                           #'O-CO_s':'initial_state',
                           #'O-O_s':'final_state',
                           }

    # set up calculation
    # surfaces = ['CoNC', 'CoNNC', 'CoNNCC', 'FeNNC', 'FeNNCC', 'FeNNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'VNCC', 'VNNC', 'VNNCC', 'VNNNCC', 'VNNNNCC', 'MnNC', 'MnNCC', 'MnNNC', 'MnNNCC', 'MnNNNCC', 'MnNNNNCC', 'RhNC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNC', 'RuNNNCC', 'RuNNNNCC', 'TiNC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNC', 'CrNCC', 'CrNNC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']
    surfaces = ['Pt', 'Pd', 'Cu', 'Ag', 'Au']

    inputs = {
        'code': catmap_code,
        'energies': energies,
        'rxn_expressions':List(list=[
                        'CO2_g + + H_g + ele_g + *_s <-> ^0.01eV_s <-> COOH_s',
                        'COOH_s + H_g + ele_g <-> CO_s + H2O_g', 
                        'CO_s <-> CO_g + *_s',
                        'H2_g <-> H2_g',
        ]), 
        'surface_names':List(list=surfaces), 
        'descriptor_names':List(list=['CO_s','COOH_s']), 
        'descriptor_ranges':List(list=[[-1.5, 1.5], [-1.5, 1.5]]), 
        'resolution':Int(30), 
        'voltage':Float(potential),
        'pH': Float(pH),
        'temperature':Float(300), 
        'species_definitions':Dict(dict=species_definitions), 
        'gas_thermo_mode':Str('ideal_gas'), 
        'adsorbate_thermo_mode':Str('harmonic_adsorbate'),
        'electrochemical_thermo_mode':List(list=['simple_electrochemical']),
        'scaling_constraint_dict':Dict(dict=scaling_constraint_dict), 
        'decimal_precision':Int(150), 
        'tolerance':Float(1e-20), 
        'max_rootfinding_iterations':Int(100), 
        'max_bisections':Int(3), 
        #'numerical_solver':Str('coverages'),
        'metadata': {
            'label':label, #'Calculation with CO2 assumed solvation and corrected dipole',
            'description': "Submissions with AiiDA CatMAP",
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    future = engine.submit(CalculationFactory('catmap'), **inputs)
    grp = Group.get(label='mkm_doped_graphene/TM_test_calcs')
    grp.add_nodes(future)
    # result = engine.run(CalculationFactory('catmap'), **inputs)

    # computed_result = result['log'].get_content()


@click.command()
@click.option('--potential')
@click.option('--ph')
@click.option('--energy_pk')
@click.option('--label')
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code, potential, ph, energy_pk, label):
    """Run example.

    Example usage: $ ./example_01.py --code diff@localhost

    Alternative (creates diff@localhost-test code): $ ./example_01.py

    Help: $ ./example_01.py --help
    """
    run_calculation(code, potential, ph, energy_pk, label)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
