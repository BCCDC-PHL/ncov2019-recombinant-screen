#!/usr/bin/env python3

import argparse
import functools

def by_library_id(a, b):
    library_id_a = a[0]
    cid_a, plate_a, index_set_id_a, well_a = library_id_a.split('-')
    library_id_b = b[0]
    cid_b, plate_b, index_set_id_b, well_b = library_id_b.split('-')
    if plate_a > plate_b:
        return 1
    if plate_a == plate_b:
        if well_a > well_b:
            return 1
    if plate_a == plate_b and well_a == well_b:
        return 0

    return -1


def main(args):
    lines = []
    header = None
    with open(args.input, 'r') as f:
        if not args.no_header:
            header = next(f).strip().split(args.delimiter)

        for line in f:
            lines.append(line.strip().split(args.delimiter))

    sorted_lines = sorted(lines, key=functools.cmp_to_key(by_library_id))

    if not args.no_header:
        print(args.delimiter.join(header))
    for line in sorted_lines:
        print(args.delimiter.join(line))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('-d', '--delimiter', default='\t')
    parser.add_argument('--no-header', action='store_true')
    args = parser.parse_args()
    main(args)
