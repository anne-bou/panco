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
# from panco.descriptor.flow import Flow
# from panco.descriptor.server import Server
# from panco.descriptor.curves import TokenBucket, RateLatency
from panco.lpSolvePath import  LPSOLVEPATH

def delta_from_deadlines(deadlines):
    n = len(deadlines)
    return [[deadlines[j] - deadlines[i] for j in range(n)] for i in range(n)]


class EdfModular:
    def __init__(self, network: Network, deadlines, filename="edf_tfa.lp", cst=10000):
        self.network = network
        self.deadlines = deadlines
        self.constant = cst
        self.delta = delta_from_deadlines(deadlines)
        self.filename = filename
        self.current_min_deadlines = np.zeros(network.num_flows)
        self.current_sigma = [self.network.flows[i].arrival_curve[0].sigma for i in range(self.network.num_flows)]

    def edf_delay(self, j, foi):
        # print('server : ',j,'flow : ', foi)
        if foi not in self.network.flows_in_server[j]:
            raise Exception("Flow of interest does not cross the server")
        # print(self.deadlines, self.current_min_deadlines)
        delta = [self.deadlines[foi] - self.current_min_deadlines[i] for i in range(self.network.num_flows)]
        delta[foi] = 0
        # print(delta)
        file = open(self.filename, 'w')
        # Objective
        file.write('max: t - u{};\n'.format(foi))
        # time constraints
        for i in self.network.flows_in_server[j]:
            file.write('u{0} <= t;\n'.format(i))
            file.write('u{0} <= u{1} + {2};\n'.format(i, foi, delta[i]))
            file.write('u{0} >= s - {1} b{0};\n'.format(i, self.constant))
            file.write('u{0} <= s + {1} - {1}b{0};\n'.format(i, self.constant))
        # edf constraints
        for i in self.network.flows_in_server[j]:
            file.write('Au{0} = Dt{0};\n'.format(i))
        # service constraint
        for i in self.network.flows_in_server[j]:
            file.write('Dt{0} - As{0} + '.format(i))
        file.write('{0} >= {1}t - {1}s;\n'.format(self.network.servers[j].service_curve[0].latency *
                                                  self.network.servers[j].service_curve[0].rate,
                                                  self.network.servers[j].service_curve[0].rate))
        for i in self.network.flows_in_server[j]:
            file.write('Dt{0} - As{0} + '.format(i))
        file.write('0 >= 0;\n')
        # arrival constraint
        for i in self.network.flows_in_server[j]:
            file.write('Au{0} - As{0} <= {1} + {2}u{0} - {2}s + {3}b{0};\n'.
                       format(i, self.current_sigma[i],
                              self.network.flows[i].arrival_curve[0].rho, self.constant))
            file.write('Au{0} - As{0} <= {1} - {1}b{0};\n'.
                       format(i, self.constant))
            file.write('Au{0} - As{0} >= 0;\n'.format(i))
        # boolean constraints
        file.write('b{0} = 0;\n'.format(foi))
        for i in self.network.flows_in_server[j]:
            file.write('b{0} <= 1;\n'.format(i))
        for i in self.network.flows_in_server[j]:
            file.write('int b{0};\n'.format(i))
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE,
                   encoding='utf-8').stdout
        # print(s)
        return float(s.split()[-1])

    def edf_backlog(self, j, foi):
        # print('server : ',j,'flow : ', foi)
        if foi not in self.network.flows_in_server[j]:
            raise Exception("Flow of interest does not cross the server")
        # print(self.deadlines, self.current_min_deadlines)
        delta = [self.deadlines[foi] - self.current_min_deadlines[i] for i in range(self.network.num_flows)]
        delta[foi] = 0
        # print(delta)
        file = open(self.filename, 'w')
        # Objective
        file.write('max: At{0} - Dt{0};\n'.format(foi))
        # time constraints
        for i in self.network.flows_in_server[j]:
            file.write('u{0} <= t;\n'.format(i))
            file.write('u{0} <= u{1} + {2};\n'.format(i, foi, delta[i]))
            file.write('u{0} >= s - {1} b{0};\n'.format(i, self.constant))
            file.write('u{0} <= s + {1} - {1}b{0};\n'.format(i, self.constant))
        # edf constraints
        for i in self.network.flows_in_server[j]:
            file.write('Au{0} = Dt{0};\n'.format(i))
        # service constraint
        for i in self.network.flows_in_server[j]:
            file.write('Dt{0} - As{0} + '.format(i))
        file.write('{0} >= {1}t - {1}s;\n'.format(self.network.servers[j].service_curve[0].latency *
                                                  self.network.servers[j].service_curve[0].rate,
                                                  self.network.servers[j].service_curve[0].rate))
        for i in self.network.flows_in_server[j]:
            file.write('Dt{0} - As{0} + '.format(i))
        file.write('0 >= 0;\n')
        # arrival constraint
        for i in self.network.flows_in_server[j]:
            file.write('Au{0} - As{0} <= {1} + {2}u{0} - {2}s + {3}b{0};\n'.
                       format(i, self.current_sigma[i],
                              self.network.flows[i].arrival_curve[0].rho, self.constant))
            file.write('Au{0} - As{0} <= {1} - {1}b{0};\n'.
                       format(i, self.constant))
            file.write('Au{0} - As{0} >= 0;\n'.format(i))
        file.write('At{0} - Au{0} <= {1} + {2}t - {2}u{0}\n;'.
                   format(foi, self.current_sigma[foi], self.network.flows[foi].arrival_curve[0].rho))
        file.write('At{0} - As{0} <= {1} + {2}t - {2}s\n;'.
                   format(foi, self.current_sigma[foi], self.network.flows[foi].arrival_curve[0].rho))
        # boolean constraints
        file.write('b{0} = 0;\n'.format(foi))
        for i in self.network.flows_in_server[j]:
            file.write('b{0} <= 1;\n'.format(i))
        for i in self.network.flows_in_server[j]:
            file.write('int b{0};\n'.format(i))
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE,
                   encoding='utf-8').stdout
        # print(s)
        return float(s.split()[-1])

    def modular_analysis(self):
        self.current_sigma = [self.network.flows[i].arrival_curve[0].sigma for i in range(self.network.num_flows)]
        delays_flows = np.zeros(self.network.num_flows)
        sigma_flows = np.zeros(self.network.num_flows)
        self.current_min_deadlines = self.deadlines
        for j in range(self.network.num_servers):
            delays_server = np.zeros(self.network.num_flows)
            for i in self.network.flows_in_server[j]:
                delays_server[i] = self.edf_delay(j, i)
                sigma_flows[i] = self.edf_backlog(j, i)
            for i in self.network.flows_in_server[j]:
                self.current_sigma[i] += delays_server[i] * self.network.flows[i].arrival_curve[0].rho
                self.current_sigma[i] = min(self.current_sigma[i], sigma_flows[i])
                delays_flows[i] += delays_server[i]
            self.current_min_deadlines = [self.current_min_deadlines[i] - delays_server[i]
                                          for i in range(self.network.num_flows)]
        return delays_flows

    def check_deadlines(self):
        delay_flows = self.modular_analysis()
        #print(delay_flows)
        for i in range(self.network.num_flows):
            if delay_flows[i] > self.deadlines[i]:
                #print(delay_flows[i], self.deadlines[i])
                return False
        return True
