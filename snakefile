configfile: "config.yaml"
#print(config)
#print(cluster_config)
from scripts.function_snakemake import parse_idxstats, check_mapping_stats, merge_bam_stats_csv


fasta_dir = config["DATA"]["ASSEMBLY"]
ref = config["DATA"]["REFERENCE"]
output_dir = config["DATA"]["OUTPUT"]
db_busco = config["TOOLS_PARAM"]["BUSCO_DATABASE"]
long_reads_dir = config["DATA"]["LONG_READS"]
short_read_dir = config["DATA"]["SHORT_READS"]
repeat_database_name = config["DATA"]["REPEAT_DB"]
script_dir = config["DATA"]["SCRIPTS"]
log_dir = f"{output_dir}LOGS/"


def get_threads(rule, default):
    """
    give threads or 'cpus-per-task from cluster_config rule : threads to SGE and cpus-per-task to SLURM
    """
    if cluster_config:
        if rule in cluster_config and 'threads' in cluster_config[rule]:
            return int(cluster_config[rule]['threads'])
        elif rule in cluster_config and 'cpus-per-task' in cluster_config[rule]:
            return int(cluster_config[rule]['cpus-per-task'])
        elif '__default__' in cluster_config and 'cpus-per-task' in cluster_config['__default__']:
            return int(cluster_config['__default__']['cpus-per-task'])
        elif '__default__' in cluster_config and 'threads' in cluster_config['__default__']:
            return int(cluster_config['__default__']['threads'])
    if workflow.global_resources["_cores"]:
        return workflow.global_resources["_cores"]
    return default

def get_fastq_file(wildcards):
    """return if file provide from cleaning or direct sample"""
    test = str(wildcards)
    split_wildcard = test.split("_")
    #print(split_wildcard)
    dico_mapping = {
                    "fasta": f"{fasta_dir}{{samples}}.fasta",
                    "index": rules.bwa_index.output
                    }
    dico_mapping.update({
                    "R1" :f"{short_read_dir}{split_wildcard[0]}_R1.fastq.gz",
                    "R2" :f"{short_read_dir}{split_wildcard[0]}_R2.fastq.gz"
                })
    # print(dico_mapping)
    return dico_mapping


def get_fastq_file_long_read(wildcards):
    """return if file provide from cleaning or direct sample"""
    test = str(wildcards)
    split_wildcard = test.split("_")
    #print(split_wildcard)
    dico_mapping = {
                    "fasta_test": f"{fasta_dir}{{samples}}.fasta"
                    }

    dico_mapping.update({
                    "reads" :f"{long_reads_dir}{split_wildcard[0]}.fastq.gz"
                })
    #print(dico_mapping)
    return dico_mapping

ASSEMBLIES, = glob_wildcards(fasta_dir+"{samples}.fasta", followlinks=True)
#print(ASSEMBLIES)
BWA_INDEX = ['amb','ann','bwt','pac','sa']


rule finale:
    input:
        report_stat_final = f"{output_dir}2_GENOME_STATS/report.html"

rule rename_contigs:
    threads: get_threads("rename_contigs",2)
    input:
        assembly = f"{fasta_dir}{{samples}}.fasta"
    output:
        sorted_fasta = f"{output_dir}1_FASTA_SORTED/{{samples}}.fasta",
    log :
        error =  f'{log_dir}fasta_sorted/fasta_sorted_{{samples}}.e',
        output = f'{log_dir}fasta_sorted/fasta_sorted_{{samples}}.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.assembly}}
                Output:
                    - Fasta_sorted: {{output.sorted_fasta}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    shell:
        f"python {script_dir}function_rename.py -i {{input.assembly}} -o {{output.sorted_fasta}} 1>{{log.output}} 2>{{log.error}}"

rule busco:
    """make quast report for our strain assembly"""
    threads: get_threads("busco", 8)
    input:
            fasta_sorted = expand(rules.rename_contigs.output.sorted_fasta, samples = ASSEMBLIES)
    output:
            busco_out = expand(f"{output_dir}2_GENOME_STATS/BUSCO/result_busco/{{samples}}.fasta/short_summary.specific.{db_busco}.{{samples}}.fasta.txt",samples = ASSEMBLIES)
    params:
            busco_out_path = f"{output_dir}2_GENOME_STATS/BUSCO/",
            busco_result = f"result_busco",
            other_option_busco = config["TOOLS_PARAM"]["BUSCO"],
            fasta = f"{output_dir}1_FASTA_SORTED/"
    log :
            error =  f'{log_dir}busco/busco.e',
            output = f'{log_dir}busco/busco.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{params.fasta}}
                Output:
                    - sa_file: {{output.busco_out}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    envmodules:
            "busco/5.1.2"
    shell:
            """
                busco -c {threads} -i {params.fasta} {params.other_option_busco} -l {db_busco} -m genome --out_path {params.busco_out_path} -o {params.busco_result} --download_path {params.busco_out_path} -f 1>{log.output} 2>{log.error}
            """

rule busco_figure:
    """make quast report for our strain assembly"""
    threads: get_threads("busco_figure", 1)
    input:
            txt_busco = rules.busco.output.busco_out
    output:
            figure = f"{output_dir}2_GENOME_STATS/BUSCO/result_busco/busco_summaries/busco_figure.png"
    params:
            folder = "busco_summaries"
    log :
            error =  f'{log_dir}busco_figure/busco.e',
            output = f'{log_dir}busco_figure/busco.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.txt_busco}}
                Output:
                    - figure_file: {{output.figure}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    envmodules:
            "busco/5.1.2"
    shell:
            """
                mkdir -p {output_dir}2_GENOME_STATS/BUSCO/result_busco/{params.folder}
                cp {input.txt_busco} {output_dir}2_GENOME_STATS/BUSCO/result_busco/{params.folder}/
                generate_plot.py -wd {output_dir}2_GENOME_STATS/BUSCO/result_busco/{params.folder}/ 1>{log.output} 2>{log.error}

            """

rule bwa_index:
    """make index with bwa for reference file"""
    threads: get_threads("bwa_index", 1)
    input:
            fasta = rules.rename_contigs.input.assembly
    output:
            expand(f'{fasta_dir}{{samples}}.fasta.{{suffix}}',samples="{samples}",suffix= BWA_INDEX)
    log :
            error =  f'{log_dir}bwa_index/{{samples}}.e',
            output = f'{log_dir}bwa_index/{{samples}}.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta strain : {{input.fasta}}
                Output:
                    - sa_file: {{output}}
                Others
                    - Threads : {{threads}}
                    - LOG error strain: {{log.error}}
                    - LOG output strain: {{log.output}}
            """
    envmodules:
            "bwa/0.7.17"
    shell:
            """
                bwa index {input.fasta} 1>{log.output} 2>{log.error}
            """
rule bwa_mem_sort_bam:
    """make bwa mem for all samples on reference"""
    threads: get_threads('bwa_mem_sort_bam', 1)
    input:
            unpack(get_fastq_file)
    output:
            bam_file =  f"{output_dir}6_MAPPING_ILLUMINA/BWA_MEM/{{samples}}.bam"
    params:
            rg = f"@RG\\tID:{{samples}}\\tSM:{{samples}}\\tPL:Illumina",
            other_options_bwa = config["TOOLS_PARAM"]["BWA_MEM"],
            other_options_samtools_view = config["TOOLS_PARAM"]["SAMTOOLS_VIEW"],
            other_options_samtools_sort = config["TOOLS_PARAM"]["SAMTOOLS_SORT"]
    log:
            error =  f'{log_dir}bwa_mem_sort_bam/{{samples}}.e',
            output = f'{log_dir}bwa_mem_sort_bam/{{samples}}.o'
    message:
            f"""
            Running {{rule}}
                Input:
                    - Fasta : {{input.fasta}}
                    - R1: {{input.R1}}
                    - R2: {{input.R2}}
                Output:
                    - Bam strain: {{output.bam_file}}
                Params:
                    - other_options_bwa: {{params.other_options_bwa}}
                    - other_options_samtools_view: {{params.other_options_samtools_view}}
                    - other_options_samtools_sort: {{params.other_options_samtools_sort}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    envmodules:
            "bwa/0.7.17",
            "samtools/1.9"
    shell:
            """
                    (bwa mem -t {threads} {input.fasta} {input.R1} {input.R2} -R '{params.rg}' |
                    samtools view -@ {threads} {params.other_options_samtools_view} |
                    samtools sort -@ {threads} {params.other_options_samtools_sort} -o {output.bam_file} ) 1>{log.output} 2>{log.error}
            """

rule samtools_index_illumina:
    """index bam for use stats"""
    threads: get_threads('samtools_index', 1)
    input:
            bam = rules.bwa_mem_sort_bam.output.bam_file
    output:
            bai = f"{output_dir}6_MAPPING_ILLUMINA/BWA_MEM/{{samples}}.bam.bai"
    log:
            error =  f'{log_dir}samtools_index/{{samples}}.e',
            output = f'{log_dir}samtools_index/{{samples}}.o'
    message:
            f"""
            Running {{rule}}
                Input:
                    - Bam : {{input.bam}}
                Output:
                    - Bai : {{output.bai}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    envmodules:
        "samtools/1.9"
    shell:
            """
                samtools index -@ {threads} {input.bam} 1>{log.output} 2>{log.error}
            """

####### Stats

rule samtools_idxstats:
    """apply samtools idxstats on all bam"""
    threads: get_threads('samtools_idxstats', 1)
    input:
            bam = rules.samtools_index_illumina.input.bam
            #bai = rules.samtools_index.output.bai
    output:
            txt_file = f"{output_dir}6_MAPPING_ILLUMINA/STATS/idxstats/{{samples}}_IDXSTATS.txt"
    log:
            error =  f'{log_dir}samtools_idxstats/{{samples}}.e',
            output = f'{log_dir}samtools_idxstats/{{samples}}.o'
    message:
            f"""
            Running {{rule}} for
                Input:
                    - Bam : {{input.bam}}
                Output:
                    - txt : {{output.txt_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    envmodules:
        "samtools/1.9"
    shell:
            """
                samtools idxstats -@ {threads} {input.bam} | tee {output.txt_file} 1>{log.output} 2>{log.error}
            """
rule merge_idxstats:
    """merge all samtools idxstats files"""
    threads : get_threads('merge_idxstats', 1)
    input :
            csv_resume = expand(rules.samtools_idxstats.output.txt_file , samples = ASSEMBLIES)

    output :
            csv_resume_merge = report(f"{output_dir}6_MAPPING_ILLUMINA/STATS/all_mapping_stats_resume.csv", category="Resume mapping info's")
    log:
            error =  f'{log_dir}merge_idxstats/all_mapping_stats_resume.e',
            output = f'{log_dir}merge_idxstats/all_mapping_stats_resume.o'
    message:
            f"""
            Running {{rule}} for
                Input:
                    - CSV_files : {{input.csv_resume}}
                Output:
                    - Merge CSV : {{output.csv_resume_merge}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    run :
        parse_idxstats(input.csv_resume, output.csv_resume_merge, sep="\t")

rule samtools_depth:
    """apply samtools depth on all bam SE end PE"""
    threads: get_threads('samtools_depth', 1)
    input:
            bam = rules.samtools_index_illumina.input.bam
            #bai = rules.samtools_index.output.bai
    output:
            txt_file = f"{output_dir}6_MAPPING_ILLUMINA/STATS/depth/{{samples}}_DEPTH.txt"
    params:
            other_options = config["TOOLS_PARAM"]["SAMTOOLS_DEPTH"]
    log:
            error =  f'{log_dir}samtools_depth/{{samples}}.e',
            output = f'{log_dir}samtools_depth/{{samples}}.o'
    message:
            f"""
            Execute {{rule}} for
                Input:
                    - Bam : {{input.bam}}
                Output:
                    - txt : {{output.txt_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    envmodules:
        "samtools/1.14"
    shell:
            """
                samtools depth -@ {threads} {params.other_options} {input.bam} | tee {output.txt_file} 1>{log.output} 2>{log.error}
            """

rule bam_stats_to_csv:
    """build csv with mean depth, median depth and mean coverage for all bam"""
    threads : get_threads('bam_stats_to_csv', 1)
    input :
            bam = rules.samtools_index_illumina.input.bam,
            bai = rules.samtools_index_illumina.output.bai,
            txt = rules.samtools_depth.output.txt_file
    output :
            csv_resume = f'{output_dir}6_MAPPING_ILLUMINA/STATS/depth/{{samples}}.csv'
    log:
            error =  f'{log_dir}bam_stats_to_csv/{{samples}}.e',
            output = f'{log_dir}bam_stats_to_csv/{{samples}}.o'
    message:
            f"""
            Running {{rule}}
                Input:
                    - Bam : {{input.bam}}
                Output:
                    - CSV : {{output.csv_resume}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    run :
        check_mapping_stats(input.bam, output.csv_resume, sep="\t")

rule merge_bam_stats:
    """merge all bam_stats_to_csv files"""
    threads : get_threads('merge_bam_stats', 1)
    input :
            csv_resume = expand(rules.bam_stats_to_csv.output.csv_resume , samples = ASSEMBLIES)
    output :
            csv_resume_merge = report(f"{output_dir}6_MAPPING_ILLUMINA/STATS/all_mapping_stats_Depth_resume.csv", category="Resume mapping info's")
    log:
            error =  f'{log_dir}merge_bam_stats/mergeResume.e',
            output = f'{log_dir}merge_bam_stats/mergeResume.o'
    message:
            f"""
            Running {{rule}}
                Input:
                    - CSV list : {{input.csv_resume}}
                Output:
                    - CSV : {{output.csv_resume_merge}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    run :
        merge_bam_stats_csv(input.csv_resume, output.csv_resume_merge, sep="\t")


rule quast_full_contigs:
    """make quast report for our strain assembly"""
    threads: get_threads("quast_full_contigs", 6)
    input:
            fasta = expand(rules.rename_contigs.output.sorted_fasta,samples=ASSEMBLIES),
            reference = ref
    output:
            quast_results = f"{output_dir}2_GENOME_STATS/QUAST/report.html"
    params:
        dir = directory(f"{output_dir}2_GENOME_STATS/QUAST/"),
        other_option_quast =config["TOOLS_PARAM"]["QUAST"]
    log :
            error =  f'{log_dir}error.e',
            output = f'{log_dir}out.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.fasta}}
                    - Ref : {{input.reference}}
                Output:
                    - sa_file: {{output.quast_results}}
                Others
                    - Threads : {{threads}}
            """
    envmodules:
            "quast/5.0.2"
    shell:
            """
                quast {input.fasta} -t {threads} {params.other_option_quast} -r {input.reference} -o {params.dir}
            """
rule minimap2:
    threads:get_threads("minimap2", 5)
    input:
        unpack(get_fastq_file_long_read),
        reference_file = ref
    output:
        bam_file = f"{output_dir}4_STRUCTURAL_VAR/minimap2/{{samples}}_sorted.bam"
    log :
        error =  f'{log_dir}minimap2/minimap2_{{samples}}.e',
        output = f'{log_dir}minimap2/minimap2_{{samples}}.o'
    params:
            option_minimap2 = config["TOOLS_PARAM"]["MINIMAP2_REF"],
            option_view = config["TOOLS_PARAM"]["SAMTOOLS_VIEW_LONG_READ_REF"],
            option_sort = config["TOOLS_PARAM"]["SAMTOOLS_SORT_LONG_READ_REF"]
    message:
        f"""
             Running {{rule}}
                Input:
                    - Reads : {{input.reads}}
                    - Reference : {{input.reference_file}}
                Output:
                    - bam_file: {{output.bam_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
        "minimap2",
        "samtools"
    shell:
        """
        minimap2 {params.option_minimap2} {input.reference_file} {input.reads} |
        samtools view -b -@ {threads} {params.option_view} |
        samtools sort -@ {threads} {params.option_sort} -o {output.bam_file}
        """

rule samtools_index:
    threads: get_threads("samtools_index",3)
    input:
        bam_file = rules.minimap2.output.bam_file
    output:
        bam_index = f"{output_dir}4_STRUCTURAL_VAR/minimap2/{{samples}}_sorted.bam.bai"
    params:
            option_index = config["TOOLS_PARAM"]["SAMTOOLS_INDEX_LONG_READ"]
    log :
        error =  f'{log_dir}minimap2/samtools_index_{{samples}}.e',
        output = f'{log_dir}minimap2/samtools_index_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - bam_file : {{input.bam_file}}
                Output:
                    - bam_index: {{output.bam_index}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
        "samtools/1.14"
    shell:
        "samtools index -@ {threads} {params.option_index} {input.bam_file} {output.bam_index}"

rule sniffles:
    threads:get_threads("sniffles",6)
    input:
        bam_file = rules.minimap2.output.bam_file,
        bam_index = rules.samtools_index.output.bam_index
    output:
        vcf_file = f"{output_dir}4_STRUCTURAL_VAR/sniffles/{{samples}}.vcf"
    params:
            option_sniffles = config["TOOLS_PARAM"]["SNIFFLES"]
    log :
        error =  f'{log_dir}sniffles/sniffles_{{samples}}.e',
        output = f'{log_dir}sniffles/sniffles_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - BAM : {{input.bam_file}}
                    - BAI : {{input.bam_index}}
                Output:
                    - vcf: {{output.vcf_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    shell:
        "sniffles -i {input.bam_file} {params.option_sniffles} -v {output.vcf_file} -t {threads}"

rule variant_calling:
    threads:get_threads("variant_calling", 1)
    input:
        vcf_file = rules.sniffles.output.vcf_file
    output:
        csv_file = f"{output_dir}4_STRUCTURAL_VAR/csv_variants/{{samples}}_variants.csv"
    log :
        error =  f'{log_dir}variants/variants_{{samples}}.e',
        output = f'{log_dir}variants/variants_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - vcf : {{input.vcf_file}}
                Output:
                    - csv_file: {{output.csv_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    shell:
        f"python {script_dir}vcf_parse.py {{input.vcf_file}} {{output.csv_file}}"

rule align_assembly:
    threads: get_threads("align_assembly",5)
    input:
        unpack(get_fastq_file_long_read),
        fasta = rules.rename_contigs.output.sorted_fasta
    output:
        bam_file = f"{output_dir}2_GENOME_STATS/COVERAGE/{{samples}}_sorted.bam"
    params:
        option_minimap2 = config["TOOLS_PARAM"]["MINIMAP2_ASSEMBLY"],
        option_view = config["TOOLS_PARAM"]["SAMTOOLS_VIEW_LONG_READ_ASSEMBLY"],
        option_sort = config["TOOLS_PARAM"]["SAMTOOLS_SORT_LONG_READ_ASSEMBLY"]
    log :
        error =  f'{log_dir}align_assembly/align_assembly_{{samples}}.e',
        output = f'{log_dir}sniffles/align_assembly_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - fasta : {{input.fasta}}
                    - reads : {{input.reads}}
                Output:
                    - bam_file: {{output.bam_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
        "minimap2",
        "samtools"
    shell:
        """
        minimap2 {params.option_minimap2} {input.fasta} {input.reads} |
        samtools view {params.option_view} -@ {threads} -b |
        samtools sort {params.option_sort} -@ {threads} -o {output.bam_file}
        """

rule coverage:
    threads: get_threads("coverage",4)
    input:
        bam_file = rules.align_assembly.output.bam_file
    output:
        coverage_file = f"{output_dir}2_GENOME_STATS/COVERAGE/{{samples}}_coverage"
    params:
        option_samtools_coverage = config["TOOLS_PARAM"]["SAMTOOLS_COVERAGE"]
    log :
        error =  f'{log_dir}coverage/coverage_{{samples}}.e',
        output = f'{log_dir}coverage/coverage_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - bam : {{input.bam_file}}
                Output:
                    - coverage: {{output.coverage_file}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
        "samtools"
    shell:
        "samtools coverage {params.option_samtools_coverage} {input.bam_file} -o {output.coverage_file}"

rule repeatmodeler:
    threads: get_threads("repeatmodeler", 12)
    input:
        fasta_files = expand(rules.rename_contigs.output.sorted_fasta, samples = ASSEMBLIES)
    output:
        database = f"{output_dir}3_REPEATMASKER/repeat_database/{repeat_database_name}-families.fa",
        fasta_concat = f"{output_dir}3_REPEATMASKER/repeat_database/concat_file.fasta"
    params:
        db_name = repeat_database_name,
        option_build_db = config["TOOLS_PARAM"]["DATABASE_REPEAT_MODELER"],
        option_repeat_modeler = config["TOOLS_PARAM"]["REPEATMODELER"],
        db_dir = f"{output_dir}3_REPEATMASKER/repeat_database/"
    log :
        error =  f'{log_dir}repeat_modeler/repeat_modeler.e',
        output = f'{log_dir}repeat_modeler/repeat_modeler.o'
    envmodules:
        "repeatmodeler/2.0.3"
    message:
        f"""
             Running {{rule}}
                Input:
                    - fasta : {{input.fasta_files}}
                Output:
                    - database: {{output.database}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    shell:
        """
        cat {input.fasta_files} > {output.fasta_concat}
        BuildDatabase -name {params.db_name} {params.option_build_db} {output.fasta_concat}
        RepeatModeler {params.option_repeat_modeler} -database {params.db_name}
        mv {params.db_name}* {params.db_dir}
        """


rule repeatmasker:
    threads: get_threads("repeatmasker", 5)
    input:
        fasta_file = rules.rename_contigs.output.sorted_fasta,
        lib_file = rules.repeatmodeler.output.database
    params:
        directory = f"{output_dir}3_REPEATMASKER/{{samples}}",
        other_option_RM = config["TOOLS_PARAM"]["REPEAT_MASKER"]
    output:
        fasta_masked = f"{output_dir}3_REPEATMASKER/{{samples}}/{{samples}}.fasta.masked"
    log :
        error =  f'{log_dir}repeatmasker/repeatmasker_{{samples}}.e',
        output = f'{log_dir}repeatmasker/repeatmasker_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.fasta_file}}
                Output:
                    - fasta_masked: {{output.fasta_masked}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    envmodules:
        "repeatmasker/4.1.2.p1"
    shell:
        "RepeatMasker -lib {input.lib_file} {params.other_option_RM} {input.fasta_file} -dir {params.directory}"


rule genome_stats:
    threads: get_threads("genome_stats",1)
    input:
        mask = rules.repeatmasker.output.fasta_masked,
        ordered = rules.rename_contigs.output.sorted_fasta,
        depth_mean = rules.coverage.output.coverage_file
    output:
        csv_stat = f"{output_dir}2_GENOME_STATS/STAT_CSV/{{samples}}.csv"
    log :
        error =  f'{log_dir}stat_csv/csv_{{samples}}.e',
        output = f'{log_dir}stat_csv/csv_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta masked : {{input.mask}}
                    - Fasta ordered : {{input.ordered}}
                    - Coverage : {{input.depth_mean}}
                Output:
                    - CSV file : {{output.csv_stat}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    shell:
        f"python {script_dir}summary_contigs.py {{input.mask}} {{input.ordered}} {{input.depth_mean}} {{output.csv_stat}}"


rule remove_contigs:
    threads: get_threads("remove_contigs",1)
    input:
        fasta_file_masked = rules.repeatmasker.output.fasta_masked,
        fasta = rules.rename_contigs.output.sorted_fasta
    output:
        final_files = f"{output_dir}5_FINAL_FASTA/{{samples}}.fasta"
    params:
        threshold = config["TOOLS_PARAM"]["REMOVE_CONTIGS_TRESHOLD"]
    log :
        error =  f'{log_dir}remove_contigs/remove_contigs_{{samples}}.e',
        output = f'{log_dir}remove_contigs/remove_contigs_{{samples}}.o'
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.fasta}}
                    - Masked_Fasta : {{input.fasta_file_masked}}
                Output:
                    - final_files: {{output.final_files}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    shell:
        f"python {script_dir}function_remove.py -m {{input.fasta_file_masked}} -f {{input.fasta}} -o {{output.final_files}} -n {{params.threshold}}"

rule mummer:
    """run mummer"""
    threads:get_threads("mummer", 3)
    input:
        final_files = rules.remove_contigs.output.final_files,
        reference_file = ref
    output:
        delta = f"{output_dir}7_ALIGNMENTS/{{samples}}.delta"

    params:
        option_nucmer = config["TOOLS_PARAM"]["NUCMER"]
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.final_files}}
                Output:
                    - deltafile: {{output.delta}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    log :
        error =  f'{log_dir}mummer/mummer_{{samples}}.e',
        output = f'{log_dir}mummer/mummer_{{samples}}.o'
    envmodules:
        "mummer4/4.0.0rc1"
    shell:
        """
            nucmer --delta {output.delta} {params.option_nucmer} {input.reference_file} {input.final_files} 1>{log.output} 2>{log.error}

        """

rule assemblytics:
    threads:get_threads("assemblytics", 1)
    input:
        delta = rules.mummer.output.delta
    params:
        prefix = f"report_{{samples}}",
        directory_as = f"{output_dir}7_ALIGNMENTS/{{samples}}_assemblytics/",
        option_assemblytics = config["TOOLS_PARAM"]["ASSEMBLYTICS"]
    output:
        dotplot = f"{output_dir}7_ALIGNMENTS/{{samples}}_assemblytics/report_{{samples}}.Assemblytics.Dotplot_filtered.png"
    message:
        f"""
             Running {{rule}}
                Input:
                    - delta : {{input.delta}}
                Output:
                    - dotplot: {{output.dotplot}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    log :
        error =  f'{log_dir}assemblytics/assemblytics_{{samples}}.e',
        output = f'{log_dir}assemblytics/assemblytics_{{samples}}.o'
    envmodules:
        "assemblytics/1.2.1"
    shell:
        """
        Assemblytics {input.delta} {params.prefix} {params.option_assemblytics}
        mv {params.prefix}* {params.directory_as}
        """

rule tapestry:
    threads: get_threads("tapestry", 8)
    input:
        unpack(get_fastq_file_long_read),
        assemblies = rules.rename_contigs.output.sorted_fasta
    output:
        tapestry_report = f"{output_dir}2_GENOME_STATS/TAPESTRY/{{samples}}/{{samples}}_tapestry_report.html"
    params:
        tapestry_dir = f"{output_dir}2_GENOME_STATS/TAPESTRY/{{samples}}/",
        other_option_tapestry = config["TOOLS_PARAM"]["TAPESTRY"]
    message:
        f"""
             Running {{rule}}
                Input:
                    - assemblies : {{input.assemblies}}
                    - long reads : {{input.reads}}
                Output:
                    - tapestry_report: {{output.tapestry_report}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}

            """
    log :
        error =  f'{log_dir}tapestry/tapestry_{{samples}}.e',
        output = f'{log_dir}tapestry/tapestry_{{samples}}.o'
    envmodules:
        "minimap2",
        "samtools"
    shell:
        """
        weave -a {input.assemblies} -r {input.reads} -c {threads} -o {params.tapestry_dir} {params.other_option_tapestry}
        mv {params.tapestry_dir}.tapestry_report.html {params.tapestry_dir}{wildcards.samples}_tapestry_report.html
        """

rule report_stats_contig:
    threads: get_threads('report_stats_contig', 1)
    input:
         csv_stat_contig = expand(rules.genome_stats.output.csv_stat, samples = ASSEMBLIES),
         figure_busco = rules.busco_figure.output.figure,
         csv_structural_variant = expand(rules.variant_calling.output.csv_file, samples = ASSEMBLIES),
         quast = rules.quast_full_contigs.output.quast_results,
         dotplot_final = expand(rules.assemblytics.output.dotplot, samples =ASSEMBLIES),
         depth_resume = rules.merge_bam_stats.output.csv_resume_merge,
         idxstats_resume = rules.merge_idxstats.output.csv_resume_merge,
         tapestry_out = expand(rules.tapestry.output.tapestry_report, samples = ASSEMBLIES)
    output:
        report = f"{output_dir}2_GENOME_STATS/report.html"
    params:
        csv_all_dir = f"{output_dir}2_GENOME_STATS/STAT_CSV",
        csv_struct_dir = f"{output_dir}4_STRUCTURAL_VAR/csv_variants",
        temporary = temp(directory(expand(f"{output_dir}2_GENOME_STATS/BUSCO/result_busco/{{samples}}.fasta/run_capnodiales_odb10", samples = ASSEMBLIES)))
    log:
            error =  f'{log_dir}report_stats/report_stats.e',
            output = f'{log_dir}report_stats/report_stats.o'
    envmodules:
        "r"
    message:
            f"""
            Running {{rule}}
                Input:
                    - csv : {{params.csv_all_dir}}
                    - figure : {{input.figure_busco}}
                Output:
                    - report : {{output.report}}
                Others
                    - Threads : {{threads}}
                    - LOG error: {{log.error}}
                    - LOG output: {{log.output}}
            """
    script:
        f"""{script_dir}report_stats_genome.Rmd"""

