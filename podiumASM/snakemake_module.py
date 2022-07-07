#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from snakemake.io import glob_wildcards
import pprint
from .global_variables import *
from .snakeWrapper import *


################################################
# GLOBAL Class

def get_list_chromosome_names(fasta_file):
    """
            Return the list of sequence name on the fasta file.
            Work with Biopython and python version >= 3.5
    """
    from Bio import SeqIO
    return [*SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))]

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


class Podium(SnakeWrapper):
    """
    to read file config
    """

    def __init__(self, workflow, config, ):
        super().__init__(workflow, config)
        # workflow is available only in __init__
        # print("\n".join(list(workflow.__dict__.keys())))
        # print(workflow.__dict__)

        # Initialisation of RattleSNP attributes
        self.fastq_path = None
        self.bam_path = None
        self.vcf_path = None
        self.samples = []
        self.run_RAXML = False
        self.run_RAXML_NG = False

        self.CHROMOSOMES = []
        self.CHROMOSOMES_WITHOUT_MITO = []

        # if provided fastq files
        self.fastq_gzip = False
        self.fastq_files_list = []
        self.fastq_files_ext = []

        self.cleaning_activated = False
        self.cleaning_tool = ""
        self.list_cleaning_tool_activated = []

        # mapping
        self.mapping_tool_activated = None
        self.mapping_activated = False
        self.mapping_stats_activated = False

        # calling
        self.calling_activated = False

        # filter
        self.vcf_filter_activated = False

       # self.__check_config_dic()

    def __build_tools_activated(self, key, allow, mandatory=False):
        tools_activate = []
        for tool, activated in self.config[key].items():
            if tool in allow:
                boolean_activated = var_2_bool(key, tool, activated)
                if boolean_activated:
                    tools_activate.append(tool)
                    self.config[key][tool] = boolean_activated
            else:
                raise ValueError(f'CONFIG FILE CHECKING FAIL for key "{key}": {tool} not avail on RattleSNP\n')
        if len(tools_activate) == 0 and mandatory:
            raise ValueError(f"CONFIG FILE CHECKING FAIL : you need to set True for at least one {key} from {allow}\n")
        return tools_activate

    def __check_fastq_files(self):
        """check if fastq file have the same extension"""
        # check if fastq file for assembly
        self.fastq_path = self.get_config_value('DATA', 'FASTQ')
        self.fastq_files_list, fastq_files_list_ext = get_files_ext(self.fastq_path, ALLOW_FASTQ_EXT)
        if not self.fastq_files_list:
            raise ValueError(
                f"CONFIG FILE CHECKING FAIL : you need to append at least on fastq with extension on "
                f"{ALLOW_FASTQ_EXT}\n")
        # check if all fastq have the same extension
        if len(fastq_files_list_ext) > 1:
            raise ValueError(
                f"CONFIG FILE CHECKING FAIL :please check 'DATA' section, key 'FASTQ', use only one extension format "
                f"from {fastq_files_list_ext} found\n")
        else:
            self.fastq_files_ext = fastq_files_list_ext[0]
        # check if fastq are gzip
        if "gz" in self.fastq_files_ext:
            self.fastq_gzip = True

    def __check_config_dic(self):
        """Configuration file checking"""
        # check output mandatory directory
        self._check_dir_or_string(level1="DATA", level2="OUTPUT")

        # check cleaning activation
        self.list_cleaning_tool_activated = self.__build_tools_activated("CLEANING", AVAIL_CLEANING)
        if len(self.list_cleaning_tool_activated) > 0:
            self.cleaning_tool = "_"+self.list_cleaning_tool_activated[0]
            self.cleaning_activated = True
        elif len(self.list_cleaning_tool_activated) > 1:
            raise ValueError(f'CONFIG FILE CHECKING FAIL for section "CLEANING": please activate only one cleaning tool avail\n')

        # check mapping activation, if not use folder name to set self.mapping_tool_activated instead of mapping tool
        self.mapping_activated = var_2_bool(tool="MAPPING", key="ACTIVATE", to_convert=self.get_config_value("MAPPING", "ACTIVATE"))
        if self.mapping_activated:
            self.mapping_stats_activated = var_2_bool(tool="MAPPING", key="BUILD_STATS", to_convert=self.get_config_value("MAPPING", "BUILD_STATS"))
            self.mapping_tool_activated = self.get_config_value("MAPPING", "TOOL")
            if self.mapping_tool_activated not in AVAIL_MAPPING:
                raise ValueError(f'CONFIG FILE CHECKING FAIL for section "MAPPING" key "TOOL": {self.mapping_tool_activated} not avail on RattleSNP\n')

        # if cleaning or mapping check fastq path and
        if self.cleaning_activated or self.mapping_activated:
            self._check_dir_or_string(level1="DATA", level2="FASTQ")
            self.__check_fastq_files()
            self.samples, = glob_wildcards(f"{self.fastq_path}{{fastq}}_R1{self.fastq_files_ext}", followlinks=True)
            for sample in self.samples:
                if not Path(f"{self.fastq_path}{sample}_R2{self.fastq_files_ext}").exists():
                    ValueError(f"DATA CHECKING FAIL : The samples '{sample}' are single-end, please only use paired data: \n")
            # check reference file
            self._check_file_or_string(level1="DATA", level2="REFERENCE_FILE")

        # check SNP calling activation:
        self.calling_activated = var_2_bool(tool="SNPCALLING", key="", to_convert=self.get_config_value(level1="SNPCALLING"))

        if not self.mapping_activated and self.calling_activated:
            self._check_dir_or_string(level1="DATA", level2="BAM")
            self.bam_path = self.get_config_value(level1="DATA", level2="BAM")
            # self.mapping_tool_activated = Path(self.bam_path).stem
            self.samples, = glob_wildcards(f"{self.bam_path}{{bam}}.bam", followlinks=True)

        # check VCF filter activation
        self.vcf_filter_activated = var_2_bool(tool="FILTER", key="", to_convert=self.get_config_value(level1="FILTER"))
        # If only VCF filtration get vcf path
        if not self.mapping_activated and not self.calling_activated and self.vcf_filter_activated:
            self._check_file_or_string(level1="DATA", level2="VCF")
            self.vcf_path = self.get_config_value(level1="DATA", level2="VCF")

        self.run_RAXML = var_2_bool(tool="RAXML", key="", to_convert=self.get_config_value(level1="RAXML"))
        self.run_RAXML_NG = var_2_bool(tool="RAXML_NG", key="", to_convert=self.get_config_value(level1="RAXML_NG"))

        # check mitochondrial name is in fasta is not Nome
        if self.cleaning_activated or self.mapping_activated or self.calling_activated:
            self.mito_name = self.get_config_value('PARAMS', 'MITOCHONDRIAL_NAME')
            self.CHROMOSOMES = get_list_chromosome_names(self.get_config_value('DATA', 'REFERENCE_FILE'))
            if self.mito_name and self.mito_name not in self.CHROMOSOMES:
                raise NameError(
                    f'CONFIG FILE CHECKING FAIL : in the "PARAMS" section, "MITOCHONDRIAL_NAME" key: the name "{self.mito_name}" is not in fasta file {self.get_config_value("DATA", "REFERENCE_FILE")}\n')
            self.CHROMOSOMES_WITHOUT_MITO = self.CHROMOSOMES.copy()
            if self.mito_name and self.mito_name in self.CHROMOSOMES:
                self.CHROMOSOMES_WITHOUT_MITO.remove(self.mito_name)

    def __repr__(self):
        return f"{self.__class__}({pprint.pprint(self.__dict__)})"

