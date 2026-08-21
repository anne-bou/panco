#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/panco

from __future__ import annotations

__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@ens.fr"
__copyright__ = "Copyright (C) 2026"
__license__ = "BSD-3"


from panco.descriptor.network import Network
from typing import Tuple, List


def times(num_servers:int, depth:List[int]) -> Tuple[List[int], List[int], int]:
    t_min = num_servers * [0]
    t_max = num_servers * [0]
    t = 0
    for j in range(num_servers - 1, -1, -1):
        t_min[j] = t + 1
        t_max[j] = t + (depth[j] + 2)
        t = t_max[j]
    return t_min, t_max, t + 1


class PLPConstraintsTree:
    # Linear analysis for fifo tree networks using only a quadratic number of time constraints
    """ Class for generating the PLP constraints for the analysis of a tree network (no cut inside the network)

    :param network: the network to analyze
    :param foi: the flow of interest, for the objective
    :param delays_flow: the delay constraints for the SFA analysis (one constraint per flow)
    :param delays_server: the delay constraints for the TFA analysis (one constraint per server)
    """
    def __init__(self, network: Network, foi: int, delays_flow=None, delays_server=None):

        self.network = network
        dates = times(self.network.num_servers, self.network.depth)
        self.t_min = dates[0]
        self.t_max = dates[1]
        self.num_dates = dates[2]
        self.foi = foi
        self.delays_flow = delays_flow
        self.delays_server = delays_server

    def time_constraints(self, f):
        f.write('\n/* Time Constraints */\n')
        f.write('t1 <= t0;\n')
        f.write('t2 <= t1;\n')
        for j in range(self.network.num_servers - 1):
            h = self.network.successors[j][0]
            tj = self.t_min[j]
            th = self.t_min[h]
            for u in range(self.network.depth[j] + 1):
                f.write('t{0} <= t{1};\n'.format(tj + u + 1, tj + u))
                f.write('t{0} <= t{1};\n'.format(tj + u, th + u))

    def arrival_constraints(self, f):
        f.write('\n/* arrival constraints */\n')
        for i in range(self.network.num_flows):
            arrival_curve = self.network.flows[i].arrival_curve
            for tb in arrival_curve:
                j = self.network.path[i][0]
                for u in range(self.t_min[j], self.t_max[j]):
                    for v in range(u + 1, self.t_max[j] + 1):
                        f.write('f{0}s{1}t{2} - f{0}s{1}t{3} <= {4} + {5} t{2} - {5} t{3};\n'.
                                format(i, j, u, v, tb.sigma, tb.rho))

    def arrival_shaping_constraints(self, f):
        f.write('\n/* arrival shaping constraints */\n')
        for k in range(len(self.network.arrival_shaping)):
            j = self.network.arrival_shaping[k][0]
            max_service = self.network.arrival_shaping[k][2]
            for i in self.network.arrival_shaping[k][1]:
                if not j == self.network.path[i][0]:
                    print('error in shaping constraints', j, self.network.path[i][0])
                    return
            for u in range(self.t_min[j], self.t_max[j]):
                for v in range(u + 1, self.t_max[j] + 1):
                    for tb in max_service:
                        f.write('0')
                        for i in self.network.arrival_shaping[k][1]:
                                f.write('+f{0}s{1}t{2} - f{0}s{1}t{3}'.format(i, j, u, v))
                        f.write('<= {0} + {1}t{2} - {1}t{3};\n'.format(tb.sigma, tb.rho, u, v))

    def monotony_constraints(self, f):
        f.write('\n/* Monotony constraints */\n')
        for i in range(self.network.num_flows):
            for j in self.network.path[i]:
                for u in range(self.t_min[j], self.t_max[j]):
                    f.write('f{0}s{1}t{2} - f{0}s{1}t{3} >= 0; \n'.format(i, j, u, u + 1))

    def fifo_constraints(self, f):
        f.write('\n/* fifo constraints */\n')
        for i in range(self.network.num_flows):
            for j in self.network.path[i]:
                if j == self.network.num_servers - 1:
                    f.write('f{0}s{1}t{3}= f{0}s{2}t{4}; \n'.format(i, j, j + 1, 1, 0))
                else:
                    h = self.network.successors[j][0]
                    for u in range(self.network.depth[j] + 1):
                        f.write('f{0}s{1}t{3} = f{0}s{2}t{4}; \n'.format(i, j, h, self.t_min[j] + u,
                                                                                 self.t_min[h] + u))

    def sfa_delay_constraints(self, f):
        f.write('\n/* SFA delay constraints */\n')
        d = self. delays_flow  # Sfa(self.network).delay()
        if d is None:
            return
        for i in range(self.network.num_flows):
            j = self.network.path[i][-1]
            h = self.network.path[i][0]
            if j == self.network.num_servers - 1:
                f.write('t{0} - t{1} <= {2};\n'.format(0, self.t_min[h], d[i]))
            else:
                j = self.network.successors[j][0]
                for k in range(self.network.depth[j] + 2):
                    f.write('t{0} - t{1} <= {2};\n'.format(self.t_min[j] + k, self.t_min[h] + k, d[i]))

    def tfa_delay_constraints(self, f):
        f.write('\n/* TFA delay constraints */\n')

        # d, s = Tfa(self.network).analysispp()
        d = self.delays_server
        if d is None:
            return
        for j in range(self.network.num_servers):
            if j == self.network.num_servers - 1:
                f.write('t{0} - t{1} <= {2};\n'.format(0, self.t_min[j], d[j]))
            else:
                h = self.network.successors[j][0]
                for k in range(self.network.depth[h] + 2):
                    f.write('t{0} - t{1} <= {2};\n'.format(self.t_min[h] + k, self.t_min[j] + k, d[j]))

    def shaping_constraints(self, f):
        f.write('\n/* Shaping constraints (e.g. maximum rate of a link)*/\n')
        for j in range(self.network.num_servers - 1):
            h = self.network.successors[j][0]
            for u in range(self.t_min[h], self.t_max[h]):
                for v in range(u + 1, self.t_max[h] + 1):
                    for tk in self.network.servers[j].max_service_curve:
                        f.write('0')
                        for i in self.network.edges[(j, h)]:  # flows_in_server[j]:
                            f.write('+ f{0}s{1}t{2} - f{0}s{1}t{3} '.format(i, h, u, v))
                        f.write('<= {3} + {0} t{1} - {0} t{2};\n'.format(tk.rho, u, v, tk.sigma))

    def service_constraints(self, f):
        f.write('\n/* Service constraints */\n')
        for j in range(self.network.num_servers):
            u = self.t_max[j]
            if j == self.network.num_servers - 1:
                v = 0
                h = self.network.num_servers
            else:
                h = self.network.successors[j][0]
                v = self.t_max[h]
            for rl in self.network.servers[j].service_curve:
                for i in self.network.flows_in_server[j]:
                    f.write('f{0}s{1}t{2} - f{0}s{3}t{4} + '.format(i, h, v, j, u))
                f.write('{0} >= {1} t{2} - {3} t{4};\n'.format(rl.rate * rl.latency, rl.rate, v, rl.rate, u))
                for i in self.network.flows_in_server[j]:
                    f.write('f{0}s{1}t{2} - f{0}s{3}t{4} + '.format(i, h, v, j, u))
                f.write('0 >= 0;\n')


    def backlog_objective(self, file):
        if True:  # self.network.path[self.foi][-1] == self.network.num_servers - 1:
            file.write(
                'max: f{0}s{1}t0 - f{0}s{2}t0;\n'.format(self.foi, self.network.flows[self.foi].path[0],
                                                             self.network.num_servers))
            j = self.network.path[self.foi][0]
            for k in range(self.t_min[j], self.t_max[j] + 1):
                file.write('f{0}s{1}t0 - f{0}s{1}t{2} <= {3} + {4}t0 - {4}t{2};\n'.
                           format(self.foi, j, k, self.network.flows[self.foi].arrival_curve[0].sigma,
                                  self.network.flows[self.foi].arrival_curve[0].rho))
        else:
            raise Exception('flow do not stop at last server\n')

    def delay_objective(self, file):
        if self.network.path[self.foi][-1] == self.network.num_servers - 1:
            file.write('max: t0 - t{};\n'.format(self.t_min[self.network.path[self.foi][0]]))
        else:
            file.write('flow do not stop at last server\n')


