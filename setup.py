#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# This file is part of the panco project.
# https://github.com/Huawei-Paris-Research-Center/Panco

__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@huawei.com"
__copyright__ = "Copyright (C) 2022, Huawei Technologies France"
__license__ = "BSD-3"

"""The setup script."""

from setuptools import setup

README = ""
try:
    with open('README.rst') as readme_file:
        README = readme_file.read()
except:
    pass

HISTORY = ""
try:
    with open('HISTORY.rst') as history_file:
        history = history_file.read()
except:
    pass

requirements = ['numpy']

setup(
    author="Anne Bouillard",
    author_email='anne.bouillard@huawei.com',
    description="Performance Analysis with Network Calculus and Optimization",
    install_requires=requirements,
    # license="BSD-3",
    long_description=README + '\n\n' + HISTORY,
    keywords='panco',
    name='panco',
    packages=['panco', 'panco.descriptor', 'panco.fifo', 'panco.staticpriorities',
              'panco.tsn', 'panco.edf', 'panco.blind'],
    url='https://github.com//anne-bou/panco',
    version='0.1.0',
    zip_safe=False,
)
