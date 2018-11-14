#! /usr/bin/env python3

from src.start_run import start_run, print_info, print_run
from src.parse_args import setup_arg_parser, read_param_table, read_batch_file
import os

def main():
    this_dir = os.path.abspath('.')

    param_table = read_param_table()
    parser = setup_arg_parser(param_table)
    args = vars(parser.parse_args())
    print_info()

    if (args['BATCH']):
        param_list = read_batch_file(args['BATCH'][0])
        for i,params in enumerate(param_list):
            os.chdir(this_dir)
            print_run(params,i)
            start_run(params)
    else:
        params = {}
        for param in param_table:
            if args[param]:
                params[param] = args[param][0]
            else:
                params[param] = param_table[param]['default']
        run = start_run(params)


if __name__ == '__main__':
    main()
