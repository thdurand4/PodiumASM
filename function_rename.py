#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Sebastien Ravel
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--fasta_input', '-i', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create tree file')
@click.option('--fasta_output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to create tree file')
def main(fasta_input, fasta_output):
    # read fasta and save to dict
    with open(fasta_input, "r") as handle:
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    # loop to extract seq len on dict[name] = len(seq)
    dico_len_seq = {}
    for seq_name in record_dict.keys():
        if seq_name not in dico_len_seq:
            dico_len_seq[seq_name] = len(record_dict[seq_name].seq)

    # loop on sorted by len to write and rename contig
    with open(fasta_output, "w") as output_handle:
        for indice, seq_name in enumerate(sorted(dico_len_seq.keys(), key=dico_len_seq.get, reverse=True), start=1):
            new_seq_name = f"contig_{indice}"
            print(seq_name, indice, new_seq_name)

            seqObj = record_dict[seq_name]
            # print(dico_seq[seq_name])
            record = SeqRecord(seqObj.seq, id=new_seq_name, name=new_seq_name, description="")
            SeqIO.write(record, output_handle, "fasta")

if __name__ == '__main__':
    main()
