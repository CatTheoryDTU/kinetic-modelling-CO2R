# -*- coding: utf-8 -*-
# pylint: skip-file
import sys
import os
import click


@click.command()
@click.option('--file')
@click.option('--label')
def main(file, label):

    SinglefileData = DataFactory('singlefile')
    energies_file = os.path.abspath(file)
    energies = SinglefileData(energies_file)
    node = energies.store_all()

    node.label = label  #'Energies at V = -0.5V SHE, energies corrected for surface charge'
    grp = Group.get(label='mkm_doped_graphene/input_energies')
    grp.add_nodes(node)

    print(f'Node of pk: {node.pk} created ')


if __name__ == '__main__':
    main()
