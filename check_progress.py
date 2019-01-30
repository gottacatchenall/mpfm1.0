#! /usr/bin/env python

import csv, os, argparse


def read_last_line(file):
    with open(file, "rb") as f:
        first = f.readline()        # Read the first line.
        f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
        while f.read(1) != b"\n":   # Until EOL is found...
            f.seek(-2, os.SEEK_CUR) # ...jump back the read byte plus one more.
        last = f.readline()
        return last

def get_folder():
    parser = argparse.ArgumentParser(description='pass the folder which contains the runs you want to aggregate')
    parser.add_argument('folder', nargs=1)
    args = parser.parse_args()
    folder = args.folder[0]
    if not os.path.exists(folder):
        print('invalid path!')
        exit(-1)
    return folder


def main():
    file = 'patches.csv'
    col = 1
    n_gen = 1000

    folder = get_folder()
    run_directories = [x[0] for x in os.walk(folder)]
    run_directories = run_directories[1:]

    for dir in run_directories:
        path = dir + '/' + file
        line = read_last_line(path)
        gen = [x.strip() for x in line.split(',')][col]
        print(dir + ': [%s / %s]' % (gen, n_gen))


if __name__ == '__main__':
    main()
