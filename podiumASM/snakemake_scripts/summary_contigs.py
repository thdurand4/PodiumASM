#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Sebastien Ravel
import click
import sys
from collections import defaultdict, OrderedDict
from Bio import SeqIO
import pandas as pd



@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--fasta_masked', '-m', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to fasta masked')
@click.option('--depth_mean', '-d', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to input fasta')
@click.option('--telomere', '-t', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to create fasta rename to contig')
@click.option('--csv_stat', '-c', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to create fasta rename to contig')
def main(fasta_masked, depth_mean, telomere, csv_stat):
    dict_compile = defaultdict(OrderedDict)

    for record in SeqIO.parse(fasta_masked, "fasta"):
        length_seq = len(record.seq)
        n_percent = (record.seq.count("N") / length_seq) * 100
        dict_compile[record.id]["Taille_en_bp"] = length_seq
        dict_compile[record.id]["Pourcentage_de_N"] = f"{n_percent:.2f}%"

    with open(depth_mean, "r") as f1:
        header = f1.readline()
        for lignes in f1:
            ligne_split = lignes.rstrip().split("\t")
            dict_compile[ligne_split[0]]["Mean_Depth"] = ligne_split[6]

    with open(telomere, "r") as f1:
        header = f1.readline()
        for lignes in f1:
            ligne_split = lignes.rstrip().split("\t")
            dict_compile[ligne_split[0]]["GC"] = ligne_split[2]
            dict_compile[ligne_split[0]]["Start_Telomere"] = ligne_split[4]
            dict_compile[ligne_split[0]]["End_Telomere"] = ligne_split[5]

    df = pd.DataFrame.from_dict(dict_compile, orient='index')
    df.reset_index(level=0, inplace=True)
    df.rename({"index": 'Contigs'}, axis='columns', inplace=True, errors="raise")
    df.to_csv(csv_stat, index=False)

    print(df)


if __name__ == '__main__':
    main()
