#! /usr/bin/env python3

import sys, os, datetime
from tabulate import tabulate
import numpy as np
import subprocess


def print_info(param_table, args):
    info = '''
    --------------------------------------------------------------------
    Metapopulation Fragmentation Model (mpfm)
    v1.0

    A individual-based model of metapopopulation dynamics in spatiotemporally stochastic environments.

    University of Colorado at Boulder
    Dept. Ecology and Evolutionary Biology | Flaxman Lab
    Michael Catchen | michael.catchen@colorado.edu

    --------------------------------------------------------------------
    '''

    print(info)

def create_run_directory(p):
    if (p):
        name = p[0]
    else:
        name = datetime.datetime.now().strftime("%Y.%m.%d_%H.%M.%S")

    path = os.path.abspath('./data/' + name)

    try:
        os.mkdir(path)
    except:
        print('A directory with that name already exists in ./data!\nPlease use a unique directory name.')
        exit(-1)
    return path

def create_ini_file(dir_path, param_table, args):
    if (dir_path[len(dir_path)-1] != '/'):
        dir_path += '/'

    ini_name = 'params.ini'
    full_init_path = dir_path + ini_name
    f = open(full_init_path, "a")


    ints = []
    doubles = []
    strings = []
    bools = []

    for arg in args:
        name = arg
        t = param_table[arg]['type']
        if (args[arg]):
            val = args[arg][0]
        else:
            val = param_table[arg]['default']

        f.write('%s %s\n'  % (name,str(val)))


def start_proc(mpfm_path):
    process = subprocess.run(mpfm_path, stdout=True, stderr=True, shell=True)
    print(process.returncode)



def start_run(param_table, args):
    print_info(param_table, args)
    mpfm_path = os.path.abspath('./bin/mpfm')

    dir_name = args['DATA_DIRECTORY'] or None
    dir_path = create_run_directory(dir_name)

    create_ini_file(dir_path, param_table, args)
    os.chdir(dir_path)
    start_proc(mpfm_path)
