# -*- coding: utf-8 -*-
# pylint: disable-all
from aiida.orm import QueryBuilder, Node, SinglefileData
from useful_classes import bcolors
from pprint import pprint
import click
from pprint import pprint


@click.command()
@click.option('--group_number')
def main(group_number):
    qb = QueryBuilder()
    qb.append(Group, filters={'id': {'==': group_number}}, tag='group')
    qb.append(SinglefileData, with_group='group')
    results = qb.all()

    for index, row in enumerate(results):
        node = row[0]
        print(
            f'{bcolors.OKBLUE} Calculation {index}: "{node.label}" {bcolors.ENDC}'
        )
        print(f'\t Node: {node.pk}')


if __name__ == '__main__':
    main()
