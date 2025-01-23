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

# from panco.descriptor.curves import TokenBucket, RateLatency
# from panco.descriptor.server import Server
from panco.descriptor.network import Network
from panco.lpSolvePath import LPSOLVEPATH
import subprocess as sp


def delta_from_deadlines(deadlines):
    n = len(deadlines)
    return [[deadlines[j] - deadlines[i] for j in range(n)] for i in range(n)]


class EdfSinkTreeLP:
    def __init__(self, network: Network, deadlines, filename="edf.lp", cst=10000):
        self.network = network
        self.deadlines = deadlines
        self.constant = cst
        self.delta = delta_from_deadlines(deadlines)
        self.filename = filename

    def time_constraints(self, foi, file):
        file.write('\n/* Time Constraints */\n')
        for j in range(self.network.num_servers):
            if self.network.is_sink(j):
                file.write('t{0} <= t{1};\n'.format(j, self.network.num_servers))
            else:
                if len(self.network.successors[j]) == 1:
                    file.write('t{0} <= t{1};\n'.format(j, self.network.successors[j][0]))
                else:
                    raise Exception('Not a tree')
        for i in range(self.network.num_flows):
            file.write('s{0} <= t{1};\n'.format(i, self.network.num_servers))
            file.write('t{0} <= s{1} + {2}b{1};\n'.format(self.network.flows[i].path[0], i, self.constant))
            file.write('s{0} <= s{1} - {2};\n'.format(i, foi, self.delta[foi][i]))

    def service_constraints(self, file):
        file.write('\n/* Service Constraints */\n')
        for j in range(self.network.num_servers):
            rl = self.network.servers[j].service_curve[0]
            if not self.network.successors[j]:
                h = self.network.num_servers
            else:
                h = self.network.successors[j][0]
            for i in self.network.flows_in_server[j]:
                file.write('f{0}s{1}t{1} - f{0}s{2}t{2} + '.format(i, h, j))
            file.write('0 >= 0;\n')
            for i in self.network.flows_in_server[j]:
                file.write('f{0}s{1}t{1} - f{0}s{2}t{2} + '.format(i, h, j))
            file.write('{0} >=  {1}t{2} - {1}t{3};\n'.format(rl.rate * rl.latency, rl.rate, h, j))

    def arrival_constraints(self, file):
        file.write('\n/* Arrival Constraints */\n')
        for i in range(self.network.num_flows):
            path = self.network.flows[i].path + [self.network.num_servers]
            tb = self.network.flows[i].arrival_curve[0]
            file.write('/* Flow {0} */\n'.format(i))
            for (k, j) in enumerate(path):
                for h in path[k + 1:]:
                    file.write('f{0}s{1}t{2} - f{0}s{1}t{3} <= {4} + {5}t{2} - {5}t{3} + {6}b{0};\n'.
                               format(i, path[0], h, j, tb.sigma, tb.rho, self.constant))
            for j in path:
                file.write('f{0}s{1}t{2} - f{0}s{1}t{1} <= {3} + {4}s{0} - {4}t{1} + {5}b{0};\n'.
                           format(i, path[0], j, tb.sigma, tb.rho, self.constant))
                file.write('f{0}s{1}t{2} - f{0}s{1}t{1} <= {3} - {3}b{0};\n'.format(i, path[0], j, self.constant))

    def deadline_constraints(self, file):
        file.write('\n/* Deadline Constraints */\n')
        for i in range(self.network.num_flows):
            file.write('f{0}s{1}t{2} = f{0}s{2}t{2};\n'.format(i, self.network.flows[i].path[0],
                                                               self.network.num_servers))

    def boolean_constraints(self, foi, file):
        file.write('\n/* Boolean Constraints */\n')
        file.write('b{0} = 0;\n'.format(foi))
        for i in range(self.network.num_flows):
            file.write('b{0} <= 1;\n'.format(i))
        for i in range(self.network.num_flows):
            file.write('int b{0};\n'.format(i))

    def delay_objective(self, foi, file):
        file.write('max: t{0} - s{1};\n'.format(self.network.num_servers, foi))

    def write_lp_edf(self, foi):
        file = open(self.filename, 'w')
        self.delay_objective(foi, file)
        self.time_constraints(foi, file)
        self.service_constraints(file)
        self.arrival_constraints(file)
        self.deadline_constraints(file)
        self.boolean_constraints(foi, file)
        file.close()

    def compute_delay(self, foi):
        self.write_lp_edf(foi)
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE,
                   encoding='utf-8').stdout
        return float(s.split()[-1])

    def check_deadlines(self):
        for i in range(self.network.num_flows):
            d = self.compute_delay(i)
            if d > self.deadlines[i]:
                return False
        return True
