#!/usr/bin/env python
import sys
import glob, os
import re
from timeit import default_timer as timer

def main():
    filename = str(sys.argv[1])
    motif_string = "[ATCG]GAA[ATCG][ATCG]TTC[ATCG][ATCG]GAA[ATCG]"
    motif = re.compile(r"[ATCG]GAA[ATCG][ATCG]TTC[ATCG][ATCG]GAA[ATCG]", flags=re.IGNORECASE)
    #revmotif = re.compile(r"[ATCG]AAG[ATCG][ATCG]CTT[ATCG][ATCG]AAG[ATCG]", flags=re.IGNORECASE)
    revmotif = re.compile(r"[ATCG]TTC[ATCG][ATCG]GAA[ATCG][ATCG]TTC[ATCG]", flags=re.IGNORECASE)

    ministart = timer()
    with open(filename, 'r') as origfasta:
        fasta = []
        appended_piece = []
        for pieces in origfasta:
            if pieces[0] == ">":
                if appended_piece == []:
                    fasta.append(pieces.strip())
                else:
                    final = "".join(appended_piece)
                    fasta.append(final)
                    appended_piece = []
                    fasta.append(pieces.strip())
            else:
                appended_piece.append(pieces.strip())
        final = "".join(appended_piece)
        fasta.append(final)
        with open("{}.motif.tsv".format(filename.rsplit(".", 1)[0]), 'w') as new_fasta:
            motif_id = "motif_id"
            hit_num = 0
            sequence_name = "sequence_name"
            start = "start"
            stop = "stop"
            strand = "strand"
            matched_sequence = "matched_sequence"
            new_fasta.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(motif_id,sequence_name,start,stop,strand,matched_sequence))
            motif_id = motif_string
            for line in fasta:
                if line[0] == ">":
                    sequence_name = line[1:].split()[0]
                else:
                    try:
                        coordinates = re.finditer(motif, line)
                        for match in coordinates:
                            start = match.span()[0]
                            stop = match.span()[1]
                            strand = "+"
                            matched_sequence = match.group()
                            new_fasta.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(motif_id,sequence_name,start,stop,strand,matched_sequence))
                            hit_num = hit_num+1
                    except:
                        print("Something wrong with matching")
                        pass
                    try:
                        revcoordinates = re.finditer(revmotif, line)
                        for revmatch in revcoordinates:
                            start = revmatch.span()[0]
                            stop = revmatch.span()[1]
                            strand = "-"
                            matched_sequence = revmatch.group()
                            new_fasta.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(motif_id,sequence_name,start,stop,strand,matched_sequence))
                            hit_num = hit_num+1
                    except:
                        print("Something wrong with reverse matching")
                        pass
            miniend = timer()
            print(filename+": ", "Number of motif hits:", hit_num, "|", miniend-ministart," seconds.")

if __name__ == '__main__':
    main()