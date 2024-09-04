#!/usr/bin/env python3

from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser(description='Convert FASTA file to GTF format.')
parser.add_argument('-i','--input', help='Input FASTA cluster file', nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin)
parser.add_argument('-o', '--out', help='Output GTF', nargs='?',
                    type=argparse.FileType('w'),
                    default=sys.stdout)

args = parser.parse_args()

for record in SeqIO.parse(args.input, "fasta"):
    print("\t".join([record.id, "fasta_to_gtf", "ORF", "1", str(len(record.seq)), ".", "+", ".", "gene_id " + record.id]), file=args.out)
