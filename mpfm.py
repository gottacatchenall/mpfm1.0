#! /usr/bin/env python3

from src.start_run import start_run
from src.parse_args import setup_arg_parser, read_param_table


def main():    
    param_table = read_param_table()
    parser = setup_arg_parser(param_table)
    args = vars(parser.parse_args())
    run = start_run(param_table, args)


if __name__ == '__main__':
    main()
