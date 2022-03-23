configfile: "config.yaml"
#print(config)
#print(cluster_config)

fasta_dir = config["DATA"]["ASSEMBLY"]
ref = config["DATA"]["REFERENCE"]
output_dir = config["DATA"]["OUTPUT"]
db_busco = config["TOOLS_PARAM"]["BUSCO_DATABASE"]
lib_dir = config["lib"]
fasta_renamed_dir = config["dir_assemblies_renamed"]
tapestry_dir = config["tapestry_dir"]
long_reads_dir = config["long_reads"]
minimap_dir = config["minimap2_dir"]
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


ASSEMBLIES, = glob_wildcards(fasta_dir+"{samples}.fasta", followlinks=True)
print(ASSEMBLIES)

rule finale:
    input:
        #dotplot_list = expand("{samples}.Assemblytics.Dotplot_filtered.png", samples = ASSEMBLIES),
        #fasta_masked_list = expand(output_dir + "{samples}.fasta.masked", samples = ASSEMBLIES),
        #final_fasta = expand("{samples}_masked.fasta", samples = ASSEMBLIES),
        fasta_renamed = expand(f"{output_dir}1_FASTA_SORTED/{{samples}}.fasta",samples=ASSEMBLIES),
        #tapestry_files = expand(tapestry_dir + "{samples}/" + "{samples}.tapestry_report.html", samples = ASSEMBLIES),
        #bam_sorted_files = expand(minimap_dir + "{samples}_sorted.bam", samples = ASSEMBLIES),
        #bam_indexes = expand(minimap_dir + "{samples}_sorted.bam.bai", samples = ASSEMBLIES),
        #vcf_list = expand(minimap_dir + "{samples}.vcf", samples = ASSEMBLIES),
	    quast_final_results = f"{output_dir}2_GENOME_STATS/QUAST/report.pdf",
	    busco_final_results = expand(f"{output_dir}2_GENOME_STATS/BUSCO/result_busco/{{samples}}.fasta/short_summary.specific.{db_busco}.{{samples}}.fasta.txt",samples=ASSEMBLIES),
        #figure_final = f"{output_dir}BUSCO/result_busco/busco_summaries/busco_figure.png"

rule rename_contigs:
    threads:1
    input:
        assembly = f"{fasta_dir}{{samples}}.fasta"
    output:
        sorted_fasta = f"{output_dir}1_FASTA_SORTED/{{samples}}.fasta"
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
        "python function_rename.py {input.assembly} {output.sorted_fasta} 1>{log.output} 2>{log.error}"

rule quast_full_contigs:
    """make quast report for our strain assembly"""
    threads: get_threads("quast_full_contigs", 6)
    input:
            fasta = expand(rules.rename_contigs.output.sorted_fasta,samples=ASSEMBLIES),
            reference = ref
    output:
            quast_results = f"{output_dir}2_GENOME_STATS/QUAST/report.pdf"
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
rule busco:
    """make quast report for our strain assembly"""
    threads: get_threads("busco", 8)
    input:
            fasta = fasta_dir
    output:
            busco_out = expand(f"{output_dir}2_GENOME_STATS/BUSCO/result_busco/{{samples}}.fasta/short_summary.specific.{db_busco}.{{samples}}.fasta.txt",samples=ASSEMBLIES)
    params:
            busco_out_path = f"{output_dir}2_GENOME_STATS/BUSCO/",
            busco_result = f"result_busco",
            other_option_busco = config["TOOLS_PARAM"]["BUSCO"]
    log :
            error =  f'{log_dir}busco/busco.e',
            output = f'{log_dir}busco/busco.o'
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.fasta}}
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
                busco -i {input.fasta} {params.other_option_busco}  -l {db_busco} -m genome --out_path {params.busco_out_path} -o {params.busco_result} --download_path {params.busco_out_path}  -f 1>{log.output} 2>{log.error}
            """
'''
rule tapestry:
    threads:3
    input:
        assemblies = "{samples}_renamed.fasta",
        reads = long_reads_dir + "{samples}.fastq.gz"
    output:
        tapestry_report = tapestry_dir + "{samples}/" + "{samples}.tapestry_report.html"
    params:
        directory = tapestry_dir + "{samples}/"
    message:
        f"""
             Running {{rule}}
                Input:
                    - assemblies : {{input.assemblies}}
                Output:
                    - tapestry_report: {{output.tapestry_report}}
                Others
                    - Threads : {{threads}}

            """
    shell:
        "weave -a {input.assemblies} -r {input.reads} -f -c 3 -o {params.directory} -t TAACCC TTAGGG"

rule mummer:
    """run mummer"""
    threads:1
    input:
        fasta_file = "{samples}_renamed.fasta",
        reference_file = ref
    params:
        prefix = "{samples}"
    output:
        delta = "{samples}.delta"
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.fasta_file}}
                Output:
                    - deltafile: {{output.delta}}
                Others
                    - Threads : {{threads}}

            """
    envmodules:
        "mummer4/4.0.0rc1"
    shell:
        """
            nucmer -p {params.prefix} {input.reference_file} {input.fasta_file}

        """

rule assemblytics:
    threads:1
    input:
        delta_files = "{samples}.delta"
    params:
        prefix = "{samples}"
    output:
        dotplot = "{samples}.Assemblytics.Dotplot_filtered.png"
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.delta_files}}
                Output:
                    - dotplot: {{output.dotplot}}
                Others
                    - Threads : {{threads}}

            """
    envmodules:
        "assemblytics/1.2.1"
    shell:
        "Assemblytics {input.delta_files} {params.prefix} 10000 50 10000"

rule minimap2:
    threads:get_threads("minimap2", 1)
    input:
        reference_file = ref,
        reads = long_reads_dir + "{samples}.fastq.gz"
    output:
        bam_file = minimap_dir + "{samples}_sorted.bam"
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.reads}}
                Output:
                    - deltafile: {{output.bam_file}}
                Others
                    - Threads : {{threads}}

            """
    envmodules:
        "minimap2",
        "samtools"
    shell:
        """
        minimap2 -ax map-ont {input.reference_file} {input.reads} |
        samtools view -b |
        samtools sort -o {output.bam_file}
        """

rule samtools_index:
    input:
        bam_file = minimap_dir + "{samples}_sorted.bam"
    output:
        bam_index = minimap_dir + "{samples}_sorted.bam.bai"
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.bam_file}}
                Output:
                    - deltafile: {{output.bam_index}}
                Others
                    - Threads : {{threads}}

            """
    envmodules:
        "samtools/1.14"
    shell:
        "samtools index {input.bam_file} {output.bam_index}"

rule sniffles:
    threads:3
    input:
        bam_file = minimap_dir + "{samples}_sorted.bam",
        bam_index = minimap_dir + "{samples}_sorted.bam.bai"
    output:
        vcf_file = minimap_dir + "{samples}.vcf"
    message:
        f"""
             Running {{rule}}
                Input:
                    - bam : {{input.bam_file}}
                Output:
                    - vcf: {{output.vcf_file}}
                Others
                    - Threads : {{threads}}

            """
    shell:
        "sniffles -i {input.bam_file} -v {output.vcf_file}"

rule repeatmasker:
    threads:1
    input:
        fasta_file = "{samples}_renamed.fasta",
        lib_file = lib_dir + "Fungi_all.fasta"
    params:
        directory = output_dir
    output:
        fasta_masked = output_dir + "{samples}.fasta.masked"
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.fasta_file}}
                Output:
                    - fasta_masked: {{output.fasta_masked}}
                Others
                    - Threads : {{threads}}

            """
    envmodules:
        "repeatmasker/4.1.2.p1"
    shell:
        "RepeatMasker -lib {input.lib_file} -no_is {input.fasta_file} -dir {params.directory}"

rule remove_contigs:
    threads:1
    input:
        fasta_file_masked = output_dir + "{samples}.fasta.masked",
        fasta = "{samples}_renamed.fasta"
    output:
        final_files = "{samples}_masked.fasta"
    params:
        threshold = 70
    message:
        f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.fasta_file_masked}}
                Output:
                    - fasta_masked: {{output.final_files}}
                Others
                    - Threads : {{threads}}

            """
    shell:
        "python function_remove.py {input.fasta_file_masked} {input.fasta} {output.final_files} {params.threshold}"

'''