configfile: "config.yaml"
#print(config)
#print(cluster_config)

fasta_dir = config["input_assemblies"]
ref = config["ref"]
output_dir = config["output_dir"]
lib_dir = config["lib"]
fasta_renamed_dir = config["dir_assemblies_renamed"]
tapestry_dir = config["tapestry_dir"]
long_reads_dir = config["long_reads"]
minimap_dir = config["minimap2_dir"]


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

rule finale:
    input:
        #dotplot_list = expand("{samples}.Assemblytics.Dotplot_filtered.png", samples = ASSEMBLIES),
        fasta_masked_list = expand(output_dir + "{samples}.fasta.masked", samples = ASSEMBLIES),
        final_fasta = expand("{samples}_masked.fasta", samples = ASSEMBLIES),
        fasta_renamed = expand("{samples}_renamed.fasta", samples = ASSEMBLIES),
        tapestry_files = expand(tapestry_dir + "{samples}/" + "{samples}.tapestry_report.html", samples = ASSEMBLIES),
        bam_sorted_files = expand(minimap_dir + "{samples}_sorted.bam", samples = ASSEMBLIES),
        bam_indexes = expand(minimap_dir + "{samples}_sorted.bam.bai", samples = ASSEMBLIES),
        vcf_list = expand(minimap_dir + "{samples}.vcf", samples = ASSEMBLIES)

rule rename_contigs:
    threads:1
    input:
        assembly = fasta_dir + "{samples}.fasta"
    output:
        file_renamed = "{samples}_renamed.fasta"
    message:
            f"""
             Running {{rule}}
                Input:
                    - Fasta : {{input.assembly}}
                Output:
                    - deltafile: {{output.file_renamed}}
                Others
                    - Threads : {{threads}}

            """
    shell:
        "python function_rename.py {input.assembly} {output.file_renamed}"

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

