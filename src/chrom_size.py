#!/usr/bin/env python

import sys

def main():
    with open(sys.argv[1], 'r') as fasta:
        with open(sys.argv[2], 'w') as new_fasta:
            line_sum = "size"
            header = "chrom"
            for line in fasta:
                if line[0] == ">":
                    new_fasta.write("{}\t{}\n".format(header,line_sum))
                    line_sum = 0
                    header = line[1:-1].split(" ")[0]
                elif line[0] == "\n":
                    new_fasta.write("{}\t{}\n".format(header,line_sum))
                else:
                    line_sum = line_sum + len(line[0:-1])
            new_fasta.write("{}\t{}\n".format(header,line_sum))

if __name__ == '__main__':
    main()