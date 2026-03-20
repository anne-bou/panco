#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/panco

from __future__ import annotations

__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@ens.fr"
__copyright__ = "Copyright (C) 2025"
__license__ = "BSD-3"

import numpy as np
import subprocess as sp
from typing import List, Tuple
from collections import defaultdict
import re
import os
from panco.descriptor.network import Network

def lex_cplex(filename):
    dict_sol = defaultdict()
    with open(filename) as f:
        for line in f:
            if 'objective ' in line:
                values = [token for token in line.split(' ') if not token == ""]
                val = [token for token in values[3].split('"')]
                delay = float(val[1])
            if 'variable ' in line:
                values = [token for token in line.split(' ') if not token == ""]
                var = [token for token in values[1].split('"')]
                val = [token for token in values[4].split('"')]
                dict_sol[var[1]] = float(val[1])
    os.remove(filename)
    return(dict_sol)

def extract_obj(filename):
    dict_sol = defaultdict()
    with open(filename) as f:
        for line in f:
            if 'objective ' in line:
                values = [token for token in line.split(' ') if not token == ""]
                val = [token for token in values[3].split('"')]
                obj = float(val[1])
            # if 'variable ' in line:
            #     values = [token for token in line.split(' ') if not token == ""]
            #     var = [token for token in values[1].split('"')]
            #     val = [token for token in values[4].split('"')]
            #     dict_sol[var[1]] = float(val[1])
    # print(filename, ' to remove')
    os.remove(filename)
    return obj


def sol_extraction(dict_sol, network:Network):
    tab_delays = np.zeros(network.num_servers)
    dict_bursts = defaultdict()
    dict_bursts_shaping = defaultdict()
    for s in dict_sol.keys():
        if s[0] == 'd':
            tab_delays[int(float(s[1:]))] = dict_sol[s]
        if s[0] == 'x':
            s1 = [int(a) for a in re.split("[xsl]", s)[1:]]
            dict_bursts[(s1[0], s1[1], s1[2])] = dict_sol[s]
        if s[0] == 'y':
            s1 = [int(a) for a in re.split("[ysl]", s)[1:]]
            dict_bursts_shaping[(s1[0], s1[1], s1[2])] = dict_sol[s]
    return tab_delays, dict_bursts, dict_bursts_shaping


def sol_extraction_burst(dict_sol, network:Network):
    tab_delays = np.zeros(network.num_servers)
    dict_bursts = defaultdict()
    dict_bursts_shaping = defaultdict()
    for s in dict_sol.keys():
        if s[0] == 'x':
            network.flows[int(float(s[1:]))].arrival_curve[0].sigma = dict_sol[s]
    return network


def sol_extraction_burst_lp(dict_sol, n_cut):
    tab_bursts = np.zeros(n_cut)
    for s in dict_sol.keys():
        if s[0] == 'x':
            tab_bursts[int(float(s[1:]))]= dict_sol[s]
    return tab_bursts

def cplex_solve(filename):
    # file_read = "tfa_new_cplex.lp"
    file_sol = "soly"
    file = open('cplex_cmd', 'w')
    file.write ('read {} \n'.format(filename))
    file.write ('optimize\n')
    file.write('write {} sol\n'.format(file_sol))
    file.write ('quit')
    file.close()
    return sp.call(["/opt/ibm/ILOG/CPLEX_Studio2211/cplex/bin/x86-64_linux/cplex",
                    "-f", "cplex_cmd"], stdout=sp.DEVNULL)
