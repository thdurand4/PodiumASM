import pandas as pd
import click
import re
import allel
from Bio import SeqIO

@click.command(context_settings={'help_option_names': ('-h', '--help'), "max_content_width": 800})
@click.option('--vcf', '-v', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input gff file')
@click.option('--output', '-o', default=None,
              type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to output file')
@click.option('--ref_fasta', '-r', default=None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
              required=True, show_default=True, help='Path to input reference genome fasta file')


def main(vcf, output, ref_fasta):

    df = allel.vcf_to_dataframe(vcf)
    list_contigs = []
    list_del = []
    list_ins = []
    list_inv = []
    list_trans = []
    list_dup = []

    for record in SeqIO.parse(ref_fasta, "fasta"):
        list_contigs.append(str(record.id))

    ID = list((df.ID))
    chrm = list((df.CHROM))

    for element in list_contigs:

        deletions = 0
        insertion = 0
        inversion = 0
        translocations = 0
        duplications = 0

        for k in range(len(chrm)):
            if chrm[k] == element:
                variant = ID[k].split(".")
                if variant[1] == 'DEL':
                    deletions += 1
                if variant[1] == 'INS':
                    insertion += 1
                if variant[1] == 'INV':
                    inversion += 1
                if variant[1] == 'BND':
                    translocations += 1
                if variant[1] == 'DUP':
                    duplications += 1
        list_del.append(deletions)
        list_ins.append(insertion)
        list_dup.append(duplications)
        list_trans.append(translocations)
        list_inv.append(inversion)

    resume = pd.DataFrame()
    resume["contigs"] = list_contigs
    resume["deletions"] = list_del
    resume["insertions"] = list_ins
    resume["inversions"] = list_inv
    resume["translocations"] = list_trans
    resume["duplications"] = list_dup

    print(resume)

    resume.to_csv(output)

if __name__ == '__main__':
    main()