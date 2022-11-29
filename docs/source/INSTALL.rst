.. contents:: Table of Contents
   :depth: 2
   :backlinks: entry

Requirements
============

PodiumASM requires |PythonVersions| and |RVersions|.

PodiumASM is developed to work on an HPC distributed cluster.

------------------------------------------------------------------------

Install PodiumASM PyPI package
===============================

First, install the PodiumASM python package with pip.

.. code-block:: bash
   
   git clone https://github.com/thdurand4/PodiumASM.git
   cd PodiumASM
   python3 -m pip install .
   podiumASM --help

Now, follow this documentation according to what you want, local or HPC mode.

------------------------------------------------------------------------

Steps for LOCAL installation
============================

Install PodiumASM in a *local* (single machine) mode using ``podiumASM install_local`` command line.

.. click:: podiumASM.main:install_local
   :prog: podiumASM install_local
   :show-nested:


------------------------------------------------------------------------

Steps for HPC distributed cluster installation
==============================================

PodiumASM uses any available snakemake profiles to ease cluster installation and resources management.
Run the command `podiumASM install_cluster` to install on a HPC cluster.
We tried to make cluster installation as easy as possible, but it is somehow necessary to adapt a few files according to your cluster environment.


.. click:: podiumASM.main:install_cluster
   :prog: podiumASM install_cluster
   :show-nested:

1. Adapt `profile` and `cluster_config.yaml`
---------------------------------------------
f
Now that PodiumASM is installed, it proposes default configuration files, but they can be modified. Please check and adapt these files to your own system architecture.

1. Adapt the pre-formatted `f –env si`snakemake profile`` to configure your cluster options.
See the section :ref:`1. Snakemake profiles` for details.

2. Adapt the :file:`cluster_config.yaml` file to manage cluster resources such as partition, memory and threads available for each job.
See the section :ref:`2. Adapting *cluster_config.yaml*` for further details.


2. Adapt `tools_path.yaml`
--------------------------

As PodiumASM uses many tools, you must install them using env modules possibilities :

1. Using the ``module load`` mode,

.. code-block:: bash

   podiumASM install_cluster --help
   podiumASM install_cluster --scheduler slurm --env modules


Adapt the file :file:``tools_path.yaml`` - in YAML (Yet Another Markup Language) - format to indicate podiumASM where the different tools are installed on your cluster.
See the section :ref:`3. How to configure tools_path.yaml` for details.


------------------------------------------------------------------------

Advance installation
====================


1. Snakemake profiles
---------------------

The Snakemake-profiles project is an open effort to create configuration profiles allowing to execute Snakemake in various computing environments
(job scheduling systems as Slurm, SGE, Grid middleware, or cloud computing), and available at https://github.com/Snakemake-Profiles/doc.

In order to run PodiumASM on HPC cluster, we take advantages of profiles.

Quickly, see `here <https://github.com/thdurand4/PodiumASM/blob/main/podiumASM/install_files/cluster_config_SLURM.yaml>`_ an example of the Snakemake SLURM profile we used for the French national bioinformatics infrastructure at IFB.

More info about profiles can be found here https://github.com/Snakemake-Profiles/slurm#quickstart.

Preparing the profile's *config.yaml* file
******************************************

Once your basic profile is created, to finalize it, modify as necessary the ``podiumASM/podiumASM/default_profile/config.yaml`` to customize Snakemake parameters that will be used internally by PodiumASM:

.. code-block:: ini

   restart-times: 0
   jobscript: "slurm-jobscript.sh"
   cluster: "slurm-submit.py"
   cluster-status: "slurm-status.py"
   max-jobs-per-second: 1
   max-status-checks-per-second: 10
   local-cores: 1
   jobs: 200                   # edit to limit the number of jobs submitted in parallel
   latency-wait: 60000000
   use-envmodules: true        # adapt True/False for env of singularuty, but only active one possibility !
   use-singularity: false      # if False, please install all R packages listed in tools_config.yaml ENVMODULE/R
   rerun-incomplete: true
   printshellcmds: true


2. Adapting *cluster_config.yaml*
----------------------------------

In the ``cluster_config.yaml`` file, you can manage HPC resources, choosing partition, memory and threads to be used by default,
or specifically, for each rule/tool depending on your HPC Job Scheduler (see `there <https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#cluster-configuration-deprecated>`_). This file generally belongs to a Snakemake profile (see above).

.. warning::
   If more memory or threads are requested, please adapt the content
   of this file before running on your cluster.


A list of PodiumASM rules names can be found in the section :ref:`Threading rules inside podiumASM`


.. warning::
   For some rules in the *cluster_config.yaml* as `rule_graph` or `run_get_versions`,
   we use by default wildcards, please don't remove it.


3. How to configure tools_path.yaml
-----------------------------------

.. note::
    About versions of tools, the user can choose themself what version of tools to use with modules.


In the ``tools_path`` file, you can find one section: ENVMODULES. In order to fill it correctly, you have 1 options:

1. Use only ENVMODULES: in this case, fill this section with the modules available on your cluster (here is an example):

.. literalinclude:: ../../podiumASM/install_files/tools_path.yaml
   :language: YAML
   :lines: 10-18

------------------------------------------------------------------------

And more ...
-------------

Threading rules inside PodiumASM
********************************

Please find here the rules names found in PodiumASM code.
It is recommended to set threads using the snakemake command when running on a single machine,
or in a cluster configuration file to manage cluster resources through the job scheduler.
This would save users a painful exploration of the snakefiles of PodiumASM.

.. code-block:: python

   rename_contigs
   busco
   busco_figure
   bwa_index
   bwa_mem_sort_bam
   samtools_index_illumina
   samtools_idxstats
   merge_idxstats
   samtools_depth
   samtools_depth_to_csv
   merge_samtools_depth_stats
   quast_full_contigs
   minimap2
   samtools_index
   sniffles
   variant_per_contig
   align_assembly
   coverage
   repeatmodeler
   repeatmasker
   remove_contigs
   mummer
   assemblytics
   tapestry
   genome_stats
   report_stats_contig
   finale



.. |PythonVersions| image:: https://img.shields.io/badge/python-3.7%2B-blue
   :target: https://www.python.org/downloads
   :alt: Python 3.7+

.. |RVersions| image:: https://img.shields.io/badge/R-%3E%3D4.0-red
   :target: https://cran.r-project.org/
   :alt: R 4.0+

.. |SnakemakeVersions| image:: https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io
   :alt: Snakemake 5.10.0+

.. |Singularity| image:: https://img.shields.io/badge/singularity-≥3.3.0-7E4C74.svg
   :target: https://sylabs.io/docs/
   :alt: Singularity 3.10.0+

.. |graphviz| image:: https://img.shields.io/badge/graphviz-%3E%3D2.40.1-green
   :target: https://graphviz.org/
   :alt: graphviz 2.40.1+
