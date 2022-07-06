#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from pathlib import Path
from setuptools import setup, find_packages
# add for remove error with pip install -e . with pyproject.toml
import site
import sys

site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

NAME = "podiumASM"
URL = "https://github.com/sravel/podiumASM"
CURRENT_PATH = Path(__file__).resolve().parent
VERSION = str(CURRENT_PATH.joinpath("podiumASM", "VERSION").open('r').readline().strip())

__doc__ = """PodiumASM is best workflow to helped selected of best assemblies between multiple assemblers tools"""


def main():
    setup(
        # Project information
        name=NAME,
        version=VERSION,
        url=URL,
        project_urls={
            "Bug Tracker": f"{URL}/issues",
            "Documentation": f"https://{NAME}.readthedocs.io/en/latest/",
            "Source Code": URL
        },
        download_url=f"{URL}/archive/{VERSION}.tar.gz",
        author="""Ravel Sebastien (CIRAD)
                    Théo DURANT (internship M2 CIRAD).
                    Simon BACHE (internship M2 CIRAD).
                """,
        author_email="sebastien.ravel@cirad.fr",
        description=__doc__.replace("\n", ""),
        long_description=CURRENT_PATH.joinpath('README.rst').open("r", encoding='utf-8').read(),
        long_description_content_type='text/x-rst',
        license='GPLv3',

        # docs compilation utils
        command_options={
            'build_sphinx': {
                'project': ('setup.py', NAME),
                'version': ('setup.py', VERSION),
                'release': ('setup.py', VERSION),
                'source_dir': ('setup.py', CURRENT_PATH.joinpath("docs", "source").as_posix()),
                'build_dir': ('setup.py', CURRENT_PATH.joinpath("docs", "build").as_posix()),
            }},

        # Package information
        packages=find_packages(),
        include_package_data=True,
        use_scm_version={
            "version_scheme": 'release-branch-semver',
            'local_scheme': "node-and-date",
            'normalize': True,
            "root": ".",
            "relative_to": __file__,
            "fallback_version": VERSION,
            "write_to": f'{NAME}/_version.py',
        },
        # setup_requires=['setuptools_scm'],
        python_requires=">=3.6",
        install_requires=[
            'PyYAML',
            'pandas',
            'matplotlib',
            'biopython',
            'numpy',
            'snakemake',
            'tqdm',
            'click>=8.0.3',
            'cookiecutter',
            'docutils < 0.18',
            'scikit-allel',
            "sniffles>=2.0.6"
        ],
        extras_require={
            'dev': ['sphinx_copybutton',
                    'sphinx_rtd_theme',
                    'sphinx_click',
                    'tox'],
        },
        entry_points={
            NAME: [f"{NAME} = __init__"],
            'console_scripts': [f"{NAME} = {NAME}.main:main",
                                f"rename_by_length = {NAME}.scripts.rename_by_length:main",
                                f"filter_TE_contig = {NAME}.scripts.filter_TE_contig:main"]},

        # Pypi information
        platforms=['unix', 'linux'],
        keywords=[
            'snakemake',
            'SNP calling',
            'workflow'
        ],
        classifiers=[
            "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
            'License :: CeCILL-C Free Software License Agreement (CECILL-C)',
            'License :: Free for non-commercial use',
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Developers',
            'Intended Audience :: End Users/Desktop',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: R',
            'Natural Language :: English',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        options={
            'bdist_wheel': {'universal': True}
        },
        zip_safe=True,  # Don't install the lib as an .egg zipfile
    )


if __name__ == '__main__':
    main()
