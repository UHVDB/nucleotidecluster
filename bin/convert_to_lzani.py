#!/usr/bin/env python

import argparse

def parse_args(args=None):
    description = "Convert kmer-db's new2all dist output to a format accetpable by LZ-ANI."
    epilog = "Example usage: python kmerdb_new2all_to_lzani.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to distance CSV created by kmer-db distance following kmer-db new2all.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output CSV reformatted to match LZ-ANI's filter file structure.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)


def modify_ani_file(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Extract the header and sequence IDs
    header = lines[0].strip().split(',')[0:-1]
    references = header[1:]  # Skip the first column (query names)

    # identify the query IDs
    hits = {}
    for line in lines[1:]:
        query_id = line.strip().split(',')[0]
        hits[query_id] = line.strip().split(query_id + ',')[1]

    # if a reference is not in hits, add it with a comma
    # else add the corresponding hit line
    mod_lines = []
    for ref in references:
        if ref not in hits:
            mod_lines.append(f"{ref},\n")
        else:
            mod_lines.append(f"{ref},{hits[ref]}\n")

    # Combine the modified header, reference rows, and original data
    final_lines = [lines[0]] + mod_lines

    # Write the modified content to the output file
    with open(output_file, 'w') as outfile:
        outfile.writelines(final_lines)


def main(args=None):
    args = parse_args(args)

    # Modify the file
    modify_ani_file(args.input, args.output)

if __name__ == "__main__":
    main()
