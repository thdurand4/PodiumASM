.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

How to create a workflow
========================

PodiumASM allows you to build a workflow using a simple ``config.yaml`` configuration file :

* First, provide the data paths
* Second, activate the requested tools for assembly and correction.
* Third, activate the tools for quality checking of assemblies.
* And last, manage the tools parameters.

To create this file, just run:


.. click:: podiumASM.main:create_config
    :prog: podiumASM create_config
    :show-nested:

Then, edit the relevant sections of the file to customize your flavor of a workflow.


1. Providing data
------------------

First, indicate the data path in the ``config.yaml`` configuration file:

.. literalinclude:: ../../podiumASM/install_files/config.yaml
    :language: YAML
    :lines: 1-8

Find here a summary table with the description of each data needed to run PodiumASM :

.. csv-table::
    :header: "Input", "Description"
    :widths: auto

    "LONG_READ", "Indicates the path to the directory with *LongRead* sequence data (fastq.gz format) to perform minimap2."
    "REFERENCE","Only one REFERENCE genome file will be used in each PodiumASM run. This REFERENCE will be used for various quality steps (i.e. ASSEMBLYTICS, QUAST)"
    "ASSEMBLY", "Provide your assembly file in one directory"
    "REPEAT_DATABASE","Provide Uniq Repeat element Database of your organism which it be used during the repeatMasker step to annotate and mask ETs in assemblies"
    "ILLUMINA", "True or False to active rules using illumina shortread"
    "SHORT_READ", "OPTIONAL : Indicates the path to the directory with *Illumina* sequence data (fastq.gz format) use paired-end data. All fastq files need to be homogeneous in their extension name. Please use *run1_R1* and *run1_R2* nomenclature."
    "OUTPUT","output *path* directory"

.. warning::

    For FASTQ, the naming conventions accepted by PodiumASM are either *NAME.fastq.gz* or *NAME.fq.gz* or *NAME.fastq* or *NAME.fq*. Use preferentially short names and avoid special characters to avoid report failure. Please do not use the long name provided directly by the sequencing machine.

    All fastq files have to be homogeneous on their extension, and can be compressed.

    Reference fasta file needs a fasta or fa extension, uncompressed.

2. Parameters for some specific tools
--------------------------------------

You can manage tools parameters on the params section in the ``config.yaml`` file.

``Busco`` specific options:

* If BUSCO is activated, you must provide to PodiumASM the path of a Busco database *OR* only the database name (See the `Busco documentation <https://busco.ezlab.org/busco_userguide.html#genome-mode-assessing-a-genome-assembly>`_).This parameter cannot be empty.

The standard parameters used in PodiumASM are shown below. Feel free to adapt it to your own requirements.

.. literalinclude:: ../../podiumASM/install_files/config.yaml
    :language: YAML
    :lines: 10-27

.. warning::
    Please check documentation of each tool (outside of PodiumASM, and make sure that the settings are correct!)


------------------------------------------------------------------------

How to run the workflow
=======================

Before attempting to run PodiumASM, please verify that you have already modified the ``config.yaml`` file as explained in :ref:`1. Providing data`.

.. warning::

   Due to a bug of CookieCutter before attempting to run PodiumASM you have to go in PodiumASM profile and comment one line in slurm-submit.py  : 

.. code-block:: bash

   cd PodiumASM/podiumASM/default_profile
   nano slurm-submit.py


If you installed PodiumASM on a HPC cluster with a job scheduler, you can run:


.. click:: podiumASM.main:run_cluster
    :prog: podiumASM run_cluster
    :show-nested:


------------------------------------------------------------------------


.. click:: podiumASM.main:run_local
    :prog: podiumASM run_local
    :show-nested:

------------------------------------------------------------------------

Advance run
===========

Providing more resources
--------------------------

If the cluster default resources are not sufficient, you can edit the ``cluster_config.yaml`` file. See :ref:`2. Adapting *cluster_config.yaml*`:

.. click:: podiumASM.main:edit_cluster_config
    :prog: podiumASM edit_cluster_config
    :show-nested:


------------------------------------------------------------------------

Output on PodiumASM
===================

The architecture of the PodiumASM output is designed as follow:

.. code-block:: bash

   OUTPUT_PODIUMASM/
   ├── 1_FASTA_SORTED
   |   ├── SAMPLE_1
   |   ├── SAMPLE_2
   |   ├── ...
   ├── 2_GENOME_STATS
   │   ├── BUSCO
   │   │   ├── file_versions.tsv
   │   │   ├── lineages
   │   │   └── result_busco     
   │   ├── COVERAGE
   |       ├── SAMPLE_1
   |       ├── SAMPLE_2
   |       ├── ...
   │   ├── QUAST
   |       ├── REPORT_QUAST
   │   ├── STAT_CSV
   |       ├── SAMPLE_1
   |       ├── SAMPLE_2
   |       ├── ...
   │   └── TAPESTRY
   |       ├── SAMPLE_1
   |       ├── SAMPLE_2
   |       ├── ...
   ├── 3_REPEATMASKER
   │       ├── SAMPLE_1
   |       ├── SAMPLE_2
   |       ├── ...
   ├── 4_STRUCTURAL_VAR
   │   ├── csv_variants
   │   ├── minimap2
   │   └── sniffles
   ├── 5_FINAL_FASTA
   │       ├── SAMPLE_1
   |       ├── SAMPLE_2
   |       ├── ...   
   ├── 6_MAPPING_ILLUMINA
   │   ├── BWA_MEM
   │   └── STATS
   ├── 7_ALIGNMENTS
   │       ├── SAMPLE_1
   |       ├── SAMPLE_2
   |       ├── ...
   ├── LOGS
   └── FINAL_REPORT


Report
======

PodiumASM generates a useful HTML report, including the versions of tools used and, for each fastq, a summary of statistics. Please have a look at |location_link| ... and enjoy !!


.. |location_link| raw:: html

    <a href="https://itrop.ird.fr/culebront_utilities/FINAL_REPORT/PodiumASM_report.html" target="_blank">example</a>


.. important::

    To visualise the report created by PodiumASM, transfer the folder ``FINAL_RESULTS`` on your local computer and open it on any web browser.
