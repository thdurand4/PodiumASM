#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Sebastien Ravel

import pandas as pd
import click
from allel import vcf_to_dataframe
from Bio import SeqIO
from collections import defaultdict, OrderedDict
from pprint import pprint as pp


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--vcf', '-v',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to input gff file')
@click.option('--output', '-o',
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to output file')
@click.option('--ref_fasta', '-r',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to input reference genome fasta file')
def main(vcf, output, ref_fasta):
    df = vcf_to_dataframe(vcf)
    print(df)

    dico_stats = defaultdict(OrderedDict)

    ID_list = list((df.ID))
    chrm_list = list((df.CHROM))

    dict_to_name = {"DEL": "deletion",
                    "INS":"insertion",
                    "INV":"inversion",
                    "BND":"translocations",
                    "DUP":"duplications"
                    }

    for chr, type_variant in zip(chrm_list, ID_list):
        variant = type_variant.split(".")[1]
        if dict_to_name[variant] not in dico_stats[chr]:
            dico_stats[chr][dict_to_name[variant]] = 0
        else:
            dico_stats[chr][dict_to_name[variant]] += 1
    dataframe_stats = pd.DataFrame.from_dict(dico_stats, orient='index')
    dataframe_stats.reset_index(level=0, inplace=True)
    dataframe_stats.rename({"index": 'Contigs'}, axis='columns', inplace=True, errors="raise")
    # print(dataframe_stats)
    with open(output, "w") as out_csv_file:
        dataframe_stats.to_csv(out_csv_file, index=False)



if __name__ == '__main__':
    main()
