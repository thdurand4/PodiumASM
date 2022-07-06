#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Sebastien Ravel
import click
from Bio import SeqIO


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--fasta_masked', '-m', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to fasta masked')
@click.option('--fasta_no_masked', '-f', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=False, help='Path to fasta without N (not masked)')
@click.option('--fasta_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create tree file')
@click.option('--max_n', '-n', type=click.FloatRange(min=0, max=1, min_open=False, max_open=False, clamp=False),
              default=0.7, show_default=True, required=False, help='remove contigs with N > value')
@click.option('--min_size', '-s', type=click.INT,
              default=10000, show_default=True, required=False, help='remove contigs with size < value')
def main(fasta_masked, fasta_no_masked, fasta_output, max_n, min_size):
    """This programme remove contig with %masked more than value `max_n`"""
    dico_record_masked = SeqIO.to_dict(SeqIO.parse(fasta_masked, "fasta"))
    dico_record = SeqIO.to_dict(SeqIO.parse(fasta_no_masked, "fasta"))
    list_keep = []
    for id_masked in dico_record_masked:
        seq_record_masked = dico_record_masked[id_masked]
        try:
            seq_record = dico_record[id_masked]
        except KeyError as e:
            raise SystemExit(f"ID {id_masked} not in file {fasta_no_masked}")
        seq = str(seq_record_masked.seq)
        n_percent = (seq.count("N") / len(seq) * 100)
        if n_percent <= max_n * 100 and len(seq_record.seq) >= min_size:
            list_keep.append(seq_record)

    SeqIO.write(list_keep, fasta_output, 'fasta')

    print(f"Remove: {len(dico_record)-len(list_keep)}\nKeep: {len(list_keep)}\nTotal: {len(dico_record)}\n")


if __name__ == '__main__':
    main()
