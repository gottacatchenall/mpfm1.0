#! /usr/bin/env python3

import argparse
import csv

def read_param_table():
    with open('param_table.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        param_table = {rows[0]:{'type':rows[1], 'default': rows[2], 'flag':rows[3], 'desc':rows[4]} for rows in csv_reader}
        return param_table

def get_description():
    description = '''

    Metapopulation Fragmentation Model (mpfm)
    v1.0

    A individual-based model of metapopopulation dynamics in spatiotemporally stochastic environments.

    University of Colorado at Boulder
    Dept. Ecology and Evolutionary Biology | Flaxman Lab
    Michael Catchen | michael.catchen@colorado.edu

    '''

    return description

def read_batch_file(file_name):
    try:
        with open(file_name) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            nruns = len(next(csv_reader))-1
            param_list = [{} for i in range(0,nruns)]

            for row in csv_reader:
                for i in range(0, nruns):
                    param_list[i][row[0]] = row[i+1]
        return param_list

    except:
        print("Failed to read batch file!")
        exit(-1)

def setup_arg_parser(param_table):
    def get_type(type):
        type_dict = {'int': int, 'double': float, 'bool': bool, 'enum': str, 'str': str}
        return type_dict[type]

    parser = argparse.ArgumentParser(description=get_description(), formatter_class=argparse.RawTextHelpFormatter)

    for param in param_table:
        name = param
        t = param_table[param]['type']
        default = param_table[param]['default']
        flag = param_table[param]['flag']
        desc = param_table[param]['desc']

        if (type != 'bool'):
            parser.add_argument(flag, nargs=1, metavar=name, dest=name, type=get_type(t), help=name + ': ' + desc)
        else:
            parser.add_argument(flag, action='store_true', help =name + ': ' + desc)

    return parser
