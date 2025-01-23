#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/panco

from __future__ import annotations

__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@huawei.com"
__copyright__ = "Copyright (C) 2022, Huawei Technologies France"
__license__ = "BSD-3"

import numpy as np
import subprocess as sp
from typing import List

from panco.descriptor.network import Network
from panco.descriptor.flow import Flow
from panco.descriptor.server import Server
from panco.descriptor.curves import TokenBucket, RateLatency
from panco.edf.edfLP  import EdfSinkTreeLP
from panco.edf.edfModular import EdfModular


# f1 = Flow([TokenBucket(1, 1)], [0, 1])
# f2 = Flow([TokenBucket(2, 2)], [0, 1])
# f3 = Flow([TokenBucket(2, 2)], [1])
# flows = [f1, f2, f3]
# sc = Server([RateLatency(10, 1)],[])
#
# net = Network([sc, sc], flows)
# deadlines = [3, 3, 3]
#
# edf_mod = EdfModular(net, deadlines)
# DELAYS_MOD = edf_mod.modular_analysis()
#
#
# edf_lp = EdfSinkTreeLP(net, deadlines)
# DELAYS_LP = [edf_lp.compute_delay(i) for i in range(net.num_flows)]
# for f in range(net.num_flows):
#     print(deadlines[f], DELAYS_MOD[f], DELAYS_LP[f])

print('** Sink-tree tandem **')

N = 5
PATHS = [[_ for _ in range(j, N)] for j in range(N)]

#print(PATHS, deadlines)

AC = [TokenBucket(3, 1) for _ in range(N)]
SC = Server([RateLatency(3, 1)], [])

FLOWS = [Flow(AC, PATHS[j]) for j in range(N)]
# SERVERS = [Server([RateLatency(4 * (j + 1), 1)], []) for j in range(N)]

#NET = Network(SERVERS, FLOWS)

N_ITER = 100

TAB_U = [0.05, 0.1, 0.15,  0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]

# TAB_U = [0.99]
TAB_U = []
for u in TAB_U:
    servers = [Server([RateLatency((j + 1) / u, 1)], []) for j in range(N)]
    net = Network(servers, FLOWS)
    n_lp = 0
    n_mod = 0
    for i in range(N_ITER):
        # deadlines = [np.random.randint(2 * (N - j) + 2,  10 * (N - j) + 2) for j in range(N)]
        deadlines = [np.random.randint(N + 1, 5 * N) for j in range(N)]
        edf_mod = EdfModular(net, deadlines)
        edf_lp = EdfSinkTreeLP(net, deadlines)
        # print([edf_lp.compute_delay(i) for i in range(N)])
        # print(edf_mod.modular_analysis())
        # print(deadlines)
        # print(edf_mod.modular_analysis(), deadlines)
        if edf_mod.check_deadlines():
            n_mod += 1
        if edf_lp.check_deadlines():
            n_lp += 1
    print(u, '\t', n_lp / N_ITER, '\t', n_mod / N_ITER)

print('*** optimal load ***')
tab_rho = []
for i in range(N_ITER):
    deadlines = [np.random.randint( N + 1, 5 * N) for j in range(N)]
    u1 = 1
    servers1 = [Server([RateLatency((j + 1) / u1, 1)], []) for j in range(N)]
    net1 = Network(servers1, FLOWS)
    b1 = EdfSinkTreeLP(net1, deadlines).check_deadlines()
    u2 = 0.0001
    servers2 = [Server([RateLatency((j + 1) / u2, 1)], []) for j in range(N)]
    net2 = Network(servers2, FLOWS)
    b2 = EdfSinkTreeLP(net2, deadlines).check_deadlines()

    if not b2:
        u_lp = 0
    else:
        while u1 - u2 > 0.01:
            um = (u1 + u2) / 2
            serversm = [Server([RateLatency((j + 1) / um, 1)], []) for j in range(N)]
            netm = Network(serversm, FLOWS)
            bm = EdfSinkTreeLP(netm, deadlines).check_deadlines()
            if bm:
                u2 = um
            else:
                u1 = um
            #print(um)
        u_lp = u2

    u1 = 1
    servers1 = [Server([RateLatency((j + 1) / u1, 1)], []) for j in range(N)]
    net1 = Network(servers1, FLOWS)
    b1 = EdfModular(net1, deadlines).check_deadlines()
    u2 = 0.0001
    servers2 = [Server([RateLatency((j + 1) / u2, 1)], []) for j in range(N)]
    net2 = Network(servers2, FLOWS)
    b2 = EdfModular(net2, deadlines).check_deadlines()

    if not b2:
        u_mod = 0.0001
    else:
        while u1 - u2 > 0.01:
            um = (u1 + u2) / 2
            serversm = [Server([RateLatency((j + 1) / um, 1)], []) for j in range(N)]
            netm = Network(serversm, FLOWS)
            bm = EdfModular(netm, deadlines).check_deadlines()
            if bm:
                u2 = um
            else:
                u1 = um
            #print(um)
        u_mod = u2

    print(u_lp, u_mod, u_lp / u_mod)
    if u_lp > 0:
        tab_rho += [u_lp / u_mod]
        print(u_lp / u_mod)