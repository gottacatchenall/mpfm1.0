#! /usr/bin/env python3

import sys, os, datetime
from tabulate import tabulate
import numpy as np
import subprocess
import time


def print_run(params, run_num, rep_num):
    print('\t--------------------------------------------------------------------')
    print('\tRUN NUMBER: %s' % run_num)
    print('\tREPLICATE NUMBER: %s' % rep_num)
    print('\tDATA_DIRECTORY: %s' % params['DATA_DIRECTORY'])


def print_info():
    info = '''
    \t--------------------------------------------------------------------
    \tMetapopulation Fragmentation Model (mpfm)
    \tv1.0

    \tA individual-based model of metapopopulation dynamics in\n\tspatiotemporally stochastic environments.

    \tUniversity of Colorado at Boulder
    \tDept. Ecology and Evolutionary Biology | Flaxman Lab
    \tMichael Catchen | michael.catchen@colorado.edu

    \t--------------------------------------------------------------------
    '''
    print(info)



def create_run_directory(p):
    path = os.path.abspath('./data/' + p)

    try:
        os.mkdir(path)
    except:
        print('\tA directory with that name already exists in ./data!\n\tPlease use a unique directory name.')
        exit(-1)
    return path

def create_ini_file(dir_path, params):
    if (dir_path[len(dir_path)-1] != '/'):
        dir_path += '/'

    ini_name = 'params.ini'
    full_init_path = dir_path + ini_name
    f = open(full_init_path, "a")

    for param in params:
        name = param
        val = params[param]
        f.write('%s,%s\n'  % (name,str(val)))


def start_proc(mpfm_path, dir_path):
    exe = (mpfm_path) + ' ' + dir_path
    process = subprocess.run(mpfm_path, stdout=True, stderr=True, shell=True)

def write_exe_name(dir_path, mpfm_path, this_dir):
    os.chdir(this_dir)

    exe = mpfm_path + ' ' + dir_path + '\n'

    with open('lb_cmd_file', 'a') as file:
        file.write(exe)


def create_batch_run_file(this_dir, params):
    mpfm_path = os.path.abspath('./bin/mpfm')
    dir_name = params['DATA_DIRECTORY'] or str(time.time())
    dir_path = create_run_directory(dir_name)
    create_ini_file(dir_path, params)
    write_exe_name(dir_path, mpfm_path, this_dir)

def start_run(params):
    mpfm_path = os.path.abspath('./bin/mpfm')

    dir_name = params['DATA_DIRECTORY']
    dir_path = create_run_directory(dir_name)

    create_ini_file(dir_path, params)
    os.chdir(dir_path)

    start_time = time.time()

    start_proc(mpfm_path, dir_path)

    elapsed_time = time.time() - start_time
    print("\n\tRUN FINISHED IN: " + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    print('\t--------------------------------------------------------------------\n')
