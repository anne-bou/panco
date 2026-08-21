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
from panco.descriptor.network import Network
from typing import List, Tuple


def times(num_servers: int, depth: List[int])->  Tuple[List[int], List[int], int]:
    t_min =  [0 for _ in range(num_servers)]
    t_max = [0 for _ in range(num_servers)]
    t = 0
    for j in range(num_servers - 1, -1, -1):
        t_min[j] = t + 1
        t_max[j] = t + 2 ** (depth[j] + 1)
        t = t_max[j]
    return t_min, t_max, t + 1


class ELPConstraintsTreeCpx:
    # Linear analysis for fifo tree networks
    """ Class for generating the ELP constraints for the analysis of a tree network (no cut inside the network)

    :param network: the network to analyze
    :param foi: the flow of interest, for the objective

    """
    def __init__(self, network: Network, foi: int):
        self.network = network
        dates = times(self.network.num_servers, self.network.depth)
        self.t_min = dates[0]
        self.t_max = dates[1]
        self.num_dates = dates[2]
        self.foi = foi

    @property
    def matrix_order(self):
        mat = np.zeros((self.num_dates, self.num_dates))
        mat[1, 0] = 1
        mat[2, 0] = 1
        mat[2, 1] = 1
        for j in range(self.network.num_servers - 2, -1, -1):
            h = self.network.successors[j][0]
            for i in range(2 ** self.network.depth[j]):
                mat[self.t_min[j] + 2 * i, self.t_min[h] + i] = 1  # fifo
                mat[self.t_min[j] + 2 * i + 1, self.t_min[h] + i] = 1  # service
                mat[self.t_min[j] + 2 * i + 1, self.t_min[j] + 2 * i] = 1  # order
                for k in range(i + 1, 2 ** self.network.depth[j]):
                    if mat[self.t_min[h] + k, self.t_min[h] + i] == 1:
                        mat[self.t_min[j] + 2 * k, self.t_min[j] + 2 * i] = 1
                        mat[self.t_min[j] + 2 * k + 1, self.t_min[j] + 2 * i + 1] = 1
                        mat[self.t_min[j] + 2 * k + 1, self.t_min[j] + 2 * i] = 1
        return mat

    def time_constraints(self, file):
        file.write('\n\\ * Time Constraints *\\ \n')
        mat = self.matrix_order
        k, n = mat.shape
        for i in range(k):
            for j in range(k):
                if mat[i, j] == 1:
                    file.write('t{0} - t{1}<= 0\n'.format(i, j))

    def arrival_constraints(self, file):
        file.write('\n\\ * arrival constraints *\\ \n')
        mat = self.matrix_order
        for i in range(self.network.num_flows):
            path = self.network.flows[i].path
            arrival_curve = self.network.flows[i].arrival_curve
            for tb in arrival_curve:
                for k in range(self.t_min[path[0]], self.t_max[path[0]]):
                    for k1 in range(k + 1, self.t_max[path[0]] + 1):
                        if mat[k1, k] == 1:
                            file.write('f{0}s{1}t{2} - f{0}s{1}t{3}  - {5} t{2} + {5} t{3}<= {4}\n'.
                                       format(i, path[0], k, k1, tb.sigma,
                                              # self.sigma[i],
                                              tb.rho))

    def monotony_constraints(self, file):
        file.write('\n\\ * Monotony constraints *\\ \n')
        mat = self.matrix_order
        for i in range(self.network.num_flows):
            for j in self.network.path[i]:
                for k in range(self.t_min[j], self.t_max[j]):
                    for h in range(k + 1, self.t_max[j] + 1):
                        if mat[h, k] == 1:
                            file.write('f{0}s{1}t{2} - f{0}s{1}t{3} >= 0 \n'.format(i, j, k, h))

    def fifo_constraints(self, file):
        file.write('\n\\ * fifo constraints *\\ \n')
        for i in range(self.network.num_flows):
            for j in self.network.path[i]:
                if j == self.network.num_servers - 1:
                    file.write('f{0}s{1}t{2} - f{3}s{4}t{5} = 0\n'.format(i, j + 1, 0,
                                                                               i, j, 1))
                else:
                    h = self.network.successors[j][0]
                    for k in range(2 ** self.network.depth[j]):
                        file.write('f{0}s{1}t{2} - f{3}s{4}t{5} = 0\n'.format(i, h, self.t_min[h] + k,
                                                                                   i, j, self.t_min[j] + 2 * k))

    def shaping_constraints(self, file):
        file.write('\n\\ * Maximum service / shaping constraints (maximum rate of the link)*\\ \n')
        for j in range(self.network.num_servers - 1):
            h = self.network.successors[j][0]
            for u in range(self.t_min[h], self.t_max[h]):
                for v in range(u + 1, self.t_max[h] + 1):
                    mat = self.matrix_order
                    if mat[v, u] == 1:
                        for tk in self.network.servers[j].max_service_curve:
                            for i in self.network.flows_in_server[j]:
                                file.write('- f{0}s{1}t{3} + f{0}s{1}t{2}'.format(i, h, u, v))
                            file.write('- {0} t{1} + {0} t{2}<= {3}\n'.format(tk.rho, u, v, tk.sigma))

    def arrival_shaping_constraints(self, f, b=False):
        f.write('\n\\ * arrival shaping constraints *\\ \n')
        mat = self.matrix_order
        for k in range(len(self.network.arrival_shaping)):
            j = self.network.arrival_shaping[k][0]
            max_service = self.network.arrival_shaping[k][2]
            for i in self.network.arrival_shaping[k][1]:
                if not j == self.network.path[i][0]:
                    print('error in shaping constraints', j, self.network.path[i][0])
                    return
            for u in range(self.t_min[j], self.t_max[j]):
                for v in range(u + 1, self.t_max[j] + 1):
                    if mat[v, u] == 1:
                        for tb in max_service:
                            for i in self.network.arrival_shaping[k][1]:
                                f.write('-f{0}s{1}t{3} + f{0}s{1}t{2}'.format(i, j, u, v))
                            f.write('<= -{1}t{2} + {1}t{3}<= {0}\n'.format(tb.sigma, tb.rho, u, v))

    def service_constraints(self, file):
        file.write('\n\\ * Service constraints *\\ \n')
        for j in range(self.network.num_servers):
            if j == self.network.num_servers - 1:
                for rl in self.network.servers[j].service_curve:
                    i = self.network.flows_in_server[j][0]
                    file.write('f{0}s{1}t{2} - f{0}s{3}t{4} '.format(i, j+1, 0, j, 2))
                    for i in self.network.flows_in_server[j][1:]:
                        file.write('+f{0}s{1}t{2} - f{3}s{4}t{5}  '.format(i, j + 1, 0,
                                                                                   i, j, 2))
                    file.write('- {1} t{2} + {3} t{4}>= -{0}\n'.format(rl.rate * rl.latency,
                                                                              rl.rate, 0,
                                                                              rl.rate, 2))
                i = self.network.flows_in_server[j][0]
                file.write('f{0}s{1}t{2} - f{0}s{3}t{4} '.format(i, j+1, 0, j, 2))
                for i in self.network.flows_in_server[j][1:]:
                    file.write('+f{0}s{1}t{2} - f{0}s{3}t{4} '.format(i, j+1, 0, j, 2))
                file.write('>= 0\n')
            else:
                for k in range(2 ** self.network.depth[j]):
                    h = self.network.successors[j][0]
                    for rl in self.network.servers[j].service_curve:
                        i = self.network.flows_in_server[j][0]
                        file.write('f{0}s{1}t{2} - f{0}s{3}t{4} '.format(i, h, self.t_min[h] + k,
                                                                                       j, self.t_min[j] + 2 * k + 1))
                        for i in self.network.flows_in_server[j][1:]:
                            file.write('+f{0}s{1}t{2} - f{0}s{3}t{4} '.format(i, h, self.t_min[h] + k,
                                                                                       j, self.t_min[j] + 2 * k + 1))
                        file.write('- {1} t{2} + {3} t{4}>= -{0}\n'.format(rl.rate * rl.latency,
                                                                                  rl.rate, self.t_min[h] + k,
                                                                                  rl.rate, self.t_min[j] + 2 * k + 1))
                    i = self.network.flows_in_server[j][0]
                    file.write('f{0}s{1}t{2} - f{0}s{3}t{4} '.format(i, h, self.t_min[h] + k,
                                                                                       j, self.t_min[j] + 2 * k + 1))
                    for i in self.network.flows_in_server[j][1:]:
                        file.write('+f{0}s{1}t{2} - f{0}s{3}t{4} '.format(i, h, self.t_min[h] + k, j,
                                                                               self.t_min[j] + 2 * k + 1))
                    file.write('>= 0\n')



    def sfa_delay_constraints(self, f):
        pass

    def tfa_delay_constraints(self, f):
        pass

    def delay_objective(self, file):
        if self.network.path[self.foi][-1] == self.network.num_servers - 1:
            file.write('maximize\n obj: ')
            file.write('t0 - t{}\n'.format(self.t_min[self.network.path[self.foi][0]]))
        else:
            file.write('flow do not stop at last server\n')
        file.write('subject to\n')
    def backlog_objective(self, file):
        if True:  # self.network.path[self.foi][-1] == self.network.num_servers - 1:
            file.write('maximize\n obj: ')
            file.write(
                'f{0}s{1}t0 - f{0}s{2}t0\n'.format(self.foi, self.network.flows[self.foi].path[0],
                                                             self.network.num_servers))
            file.write('subject to\n')
            j = self.network.path[self.foi][0]
            for k in range(self.t_min[j], self.t_max[j] + 1):
                file.write('f{0}s{1}t0 - f{0}s{1}t{2} - {4}t0 + {4}t{2}<= {3} \n'.
                           format(self.foi, j, k, self.network.flows[self.foi].arrival_curve[0].sigma,
                                  self.network.flows[self.foi].arrival_curve[0].rho))
        else:
            raise Exception('flow do not stop at last server\n')
