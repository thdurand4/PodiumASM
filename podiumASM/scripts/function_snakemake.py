from pathlib import Path
from collections import defaultdict, OrderedDict

def get_genome_size(fasta_file):
    """
            Return the list of sequence name on the fasta file.
            Work with Biopython and python version >= 3.5
    """
    from Bio import SeqIO
    return sum([len(seq) for seq in SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")).values()])


def parse_idxstats(files_list=None, out_csv=None, sep="\t"):
    from pathlib import Path
    from collections import defaultdict, OrderedDict
    import pandas as pd
    dico_mapping_stats = defaultdict(OrderedDict)
    for csv_file in files_list:
        sample = Path(csv_file).stem.split("_")[0]
        df = pd.read_csv(csv_file, sep="\t", header=None, names=["chr", "chr_size", "map_paired", "map_single"], index_col=False)
        # print(df)
        unmap = df[df.chr == '*'].map_single.values[0]
        df = df[df.chr != '*']
        map_total = df["map_single"].sum()+df["map_paired"].sum()
        size_lib = map_total+unmap
        percent = map_total/size_lib
        dico_mapping_stats[f"{sample}"]["size_lib"] = size_lib
        dico_mapping_stats[f"{sample}"]["map_total"] = map_total
        dico_mapping_stats[f"{sample}"]["percent"] = f"{percent*100:.2f}%"
    dataframe_mapping_stats = pd.DataFrame.from_dict(dico_mapping_stats, orient='index')
    dataframe_mapping_stats.reset_index(level=0, inplace=True)
    dataframe_mapping_stats.rename({"index": 'Samples'}, axis='columns', inplace=True, errors="raise")
    with open(out_csv, "w") as out_csv_file:
        # print(f"Library size:\n{dataframe_mapping_stats}\n")
        dataframe_mapping_stats.to_csv(out_csv_file, index=False, sep=sep)


def check_mapping_stats(ref, depth_file, out_csv, sep="\t"):
    from numpy import median, mean
    import pandas as pd
    genome_size = get_genome_size(ref)
    dicoResume = defaultdict(OrderedDict)

    sample = Path(depth_file).stem.split("_DEPTH")[0]
    listMap = []
    with open(depth_file, "r") as depth_file_open:
        for line in depth_file_open:
            chr, pos, depth = line.rstrip().split("\t")
            listMap.append(int(depth))
    dicoResume[sample]["Mean Depth"] = f"{mean(listMap):.2f}"
    dicoResume[sample]["Median Depth"] = f"{median(listMap):.2f}"
    dicoResume[sample]["Max Depth"] = f"{max(listMap):.2f}"
    dicoResume[sample]["Mean Genome coverage"] = f"{(len(listMap)/genome_size)*100:.2f}%"

    dataframe_mapping_stats = pd.DataFrame.from_dict(dicoResume, orient='index')
    dataframe_mapping_stats.reset_index(level=0, inplace=True)
    dataframe_mapping_stats.rename({"index": 'Samples'}, axis='columns', inplace=True, errors="raise")
    with open(out_csv, "w") as out_csv_file:
        # print(f"Library size:\n{dataframe_mapping_stats}\n")
        dataframe_mapping_stats.to_csv(out_csv_file, index=False, sep=sep)


def merge_samtools_depth_csv(csv_files, csv_file, sep="\t"):
    # dir = Path(csv_files)
    import pandas as pd
    df = (pd.read_csv(f, sep=sep) for f in csv_files)
    df = pd.concat(df)
    df.rename(columns={'Unnamed: 0': 'Samples'}, inplace=True)
    with open(csv_file, "w") as out_csv_file:
        # print(f"All CSV infos:\n{df}\n")
        df.to_csv(out_csv_file, index=False, sep=sep)


