#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .snakeWrapper import *
from .global_variables import *
from .snakemake_module import Podium, parse_idxstats, check_mapping_stats, merge_samtools_depth_csv

logo = INSTALL_PATH.joinpath('PodiumASM_logo.png').as_posix()

__version__ = "1.0.0"

__doc__ = """Long-read sequencing is a highly accurate approach that can be used to challenging genomes, such as those containing stretches of highly repetitive elements and lot of structural variant. Long read sequencing can also be used to generate de novo assembly and genome finishing applications
Lot of tools are used to make genome assembly with long reads every day and sometimes you don't know wich Assembler tool is the best for your organism.
PodiumASM is here for you ! PodiumASM is is an open-source, scalable, modulable and traceable snakemake pipeline, able to compare multiple long read assemblies obtained from multiple assemblers tools. The workflow PodiumASM can help you to choose the best assemblies among all possibilities."""

description_tools = f"""
    Welcome to PodiumASM version: {__version__}! Created on January 2022
    @author: Sebastien Ravel (CIRAD), Theo Durand (CIRAD), Simon Bache (CIRAD)
    @email: theo.durand@cirad.fr
       
    Please cite our github: {GIT_URL}
    Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
    and GPLv3 Intellectual property belongs to IRD, CIRAD and authors.
    Documentation avail at: {DOCS}
    {get_last_version(url=GIT_URL, current_version=__version__)}"""


MODULE_FILE = f"""#%Module1.0
##
## Required internal variables
set     prefix       {INSTALL_PATH.as_posix().strip()}
set     version      {__version__.strip()}
# check if install directory exist
if {{![file exists $prefix]}} {{
    puts stderr "\t[module-info name] Load Error: $prefix does not exist"
    break
    exit 1
}}
## List conflicting modules here
conflict podium
## List prerequisite modules here
module load graphviz/2.40.1
set		fullname	CulebrONT-{__version__.strip()}
set		externalurl	"\n\t{DOCS.strip()}\n"
set		description	"\n\t{__doc__.strip()}
## Required for "module help ..."
proc ModulesHelp {{}} {{
  global description externalurl
  puts stderr "Description - $description"
  puts stderr "More Docs   - $externalurl"
}}
## Required for "module display ..." and SWDB
module-whatis   "loads the [module-info name] environment"
## Software-specific settings exported to user environment
prepend-path PATH $prefix
"""


description_tools = f"""
    Welcome to PodiumASM version: {__version__}! Created on January 2022
    @author: Sebastien Ravel (CIRAD), Theo Durand, Simon Bache
    @email: sebastien.ravel@cirad.fr
       
    Please cite our github: {GIT_URL}
    Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
    and GPLv3 Intellectual property belongs to CIRAD and authors.
    Documentation avail at: {DOCS}
    {get_last_version(url=GIT_URL, current_version=__version__)}"""


MODULE_FILE = f"""#%Module1.0
##

## Required internal variables
set     prefix       {INSTALL_PATH.as_posix().strip()}
set     version      {__version__.strip()}

# check if install directory exist
if {{![file exists $prefix]}} {{
    puts stderr "\t[module-info name] Load Error: $prefix does not exist"
    break
    exit 1
}}

## List conflicting modules here
conflict podium

## List prerequisite modules here
module load graphviz/2.40.1

set		fullname	CulebrONT-{__version__.strip()}
set		externalurl	"\n\t{DOCS.strip()}\n"
set		description	"\n\t{__doc__.strip()}

## Required for "module help ..."
proc ModulesHelp {{}} {{
  global description externalurl
  puts stderr "Description - $description"
  puts stderr "More Docs   - $externalurl"
}}

## Required for "module display ..." and SWDB
module-whatis   "loads the [module-info name] environment"

## Software-specific settings exported to user environment

prepend-path PATH $prefix

"""
