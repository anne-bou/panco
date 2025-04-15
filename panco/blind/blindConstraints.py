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
from typing import List
from panco.descriptor.network import Network


class BlindConstraints:
    """
    The class to write the linear program for computing worst-case performance bounds of a flow in a network
    under blind multiplexing. The linear constraints written are related to one flow (the flow of interest), and the
    relation to other flows is made to do the complete analysis in blindLP.py.

    :param network: the network to analyze (must be a root-directed tree)
    :param foi: the flow of interest for the performance bound (will be used in the objective and some additional \\
        linear constraints)
    :param next_foi: The successor of the foi (used for the analysis of networks  that are not trees,\\
        to be able to link different pieces of the linear program
    :param list_flows: if the tree is obtained by transformation, keeps track of the original numbering of the flows\\
        in order to write the constraints properly
    """
    def __init__(self, network: Network, foi: int, next_foi=None, list_flows=None):
        """
        Constructor for the class BlindConstraints
        """
        self.network = network
        self.foi = foi
        if next_foi is None:
            self.next_foi = 0
        else:
            self.next_foi = next_foi
        if list_flows is None:
            self.list_flows = range(self.network.num_flows)
            self.is_cyclic = False
        else:
            self.list_flows = list_flows
            self.is_cyclic = True

    def time_constraints(self, f):
        """
        Writes the time constraints for the network.

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Time Constraints */\n')
        e = self.next_foi
        f.write('t{0}e{2} <= t{1}e{2};\n'.format(self.network.num_servers - 1, self.network.num_servers, e))
        for j in range(self.network.num_servers - 1):
            h = self.network.successors[j][0]
            f.write('t{0}e{2} <= t{1}e{2};\n'.format(j, h, e))

    def arrival_constraints(self, f):
        """
        Writes the arrival constraints for each flow of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* arrival constraints */\n')
        e = self.next_foi
        for i in range(self.network.num_flows):
            j = self.network.path[i][0]
            for tb in self.network.flows[i].arrival_curve:
                for (n, h1) in enumerate(self.network.path[i]):
                    if not self.network.successors[self.network.path[i][-1]]:
                        sub_p = self.network.path[i][n + 1:] + [self.network.num_servers]
                    else:
                        sub_p = self.network.path[i][n + 1:] + [self.network.successors[self.network.path[i][-1]][0]]
                    for h2 in sub_p:
                        f.write('f{0}s{1}t{2}e{6} - f{0}s{1}t{3}e{6} <= x{4} + {5} t{2}e{6} - {5} t{3}e{6};\n'.
                                format(i, j, h2, h1, self.list_flows[i], tb.rho, e))

    def arrival_shaping_constraints(self, f):
        """
        Writes the arrival shaping constraints for each flow of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* arrival shaping constraints */\n')
        e = self.next_foi
        for k in range(len(self.network.arrival_shaping)):
            j = self.network.arrival_shaping[k][0]
            max_service = self.network.arrival_shaping[k][2]
            p_max = -1
            l_max = 0
            for i in self.network.arrival_shaping[k][1]:
                if not j == self.network.path[i][0]:
                    print('error in shaping constraints', j, self.network.path[i][0])
                    return
                if len(self.network.path[i]) > l_max:
                    p_max = self.network.path[i]
                    l_max = len(p_max)
            for (n, h1) in enumerate(p_max):
                if not self.network.successors[p_max[-1]]:
                    sub_p = p_max[n + 1:] + [self.network.num_servers]
                else:
                    sub_p = p_max[n + 1:] + [self.network.successors[p_max[-1]][0]]
                for h2 in sub_p:
                    for tb in max_service:
                        f.write('0 ')
                        for i in self.network.arrival_shaping[k][1]:
                            if not self.is_cyclic or not i == self.foi:
                                f.write('+ f{0}s{1}t{2}e{4} - f{0}s{1}t{3}e{4} '.format(i, j, h2, h1, e))
                        f.write('<= {0} + {1}t{2}e{4} - {1}t{3}e{4};\n'.format(tb.sigma, tb.rho, h2, h1, e))

    def monotony_constraints(self, f):
        """
        Writes the monotony constraints for each arrival process of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Monotony constraints */\n')
        e = self.next_foi
        for i in range(self.network.num_flows):
            j = self.network.path[i][0]
            h = self.network.path[i][-1]
            if h == self.network.num_servers - 1:
                s = self.network.num_servers
            else:
                s = self.network.successors[h][0]
            path = self.network.path[i] + [s]
            for (n, h) in enumerate(self.network.path[i]):
                f.write('f{0}s{1}t{2}e{4} - f{0}s{1}t{3}e{4} <= 0; \n'.format(i, j, h, path[n + 1], e))

    def causality_constraints(self, f):
        """
        Writes the causality constraints for each arrival process of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Causality constraints */\n')
        e = self.next_foi
        for i in range(self.network.num_flows):
            j = self.network.path[i][0]
            h = self.network.path[i][-1]
            if h == self.network.num_servers - 1:
                s = self.network.num_servers
            else:
                s = self.network.successors[h][0]
            path = self.network.path[i][1:] + [s]
            for h in path:
                f.write('f{0}s{1}t{3}e{4} - f{0}s{2}t{3}e{4} >= 0; \n'.format(i, j, h, h, e))

    def service_constraints(self, f):
        """
        Writes the service constraints for each server of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Service constraints */\n')
        e = self.next_foi
        for j in range(self.network.num_servers):
            if j == self.network.num_servers - 1:
                h = self.network.num_servers
            else:
                h = self.network.successors[j][0]
            for rl in self.network.servers[j].service_curve:
                for i in self.network.flows_in_server[j]:
                    f.write('f{0}s{1}t{1}e{3} - f{0}s{2}t{2}e{3} + '.format(i, h, j, e))
                f.write('{0} >= {1} t{2}e{4} - {1} t{3}e{4};\n'.format(rl.rate * rl.latency, rl.rate, h, j, e))
                for i in self.network.flows_in_server[j]:
                    f.write('f{0}s{1}t{1}e{3} >= f{0}s{2}t{2}e{3};\n '.format(i, h, j, e))

    def greedy_shaping_constraints(self, f):
        """
        Writes the greedy-shaping constraints for each server of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Greedy Shaping constraints */\n')
        e = self.next_foi
        for j in range(self.network.num_servers - 1):
            h = self.network.successors[j][0]
            for tk in self.network.servers[j].max_service_curve:
                f.write('0')
                for i in self.network.edges[(j, h)]:  # flows_in_server[j]:
                    if not self.is_cyclic or not i == self.foi:
                        f.write('+ f{0}s{1}t{1}e{3} - f{0}s{2}t{2}e{3}'.format(i, h, j, e))
                f.write('<= {4} + {0} t{1}e{3} - {0} t{2}e{3};\n'.format(tk.rho, h, j, e, tk.sigma))

    def fix_point_constraints(self, file):
        """
        Writes the constraints to make relation between the different pieces of the network (decomposed into a forest).
        Basically, it relates the bursts parameters of the different pieces.

        :param file: file where constraints are written
        :return: nothing
        """
        file.write('\n/* the x burst constraints*/\n')
        file.write('x{4} <= f{0}s{1}t{2}e{3} - f{0}s{2}t{2}e{3};\n'.
                   format(self.foi, self.network.path[self.foi][0],
                          self.network.path[self.foi][-1] + 1, self.next_foi, self.next_foi))

    def _paths_from_server(self, server, group_flows):
        max_server = 0
        for i in group_flows:
            if not self.network.flows[i].path[0] == server:
                raise Exception("groups of flows do not all start at the right server")
            max_server = max(max_server, self.network.flows[i].path[-1])
        list_paths = [[server]]
        server_cur = server
        while server_cur < max_server:
            h = self.network.successors[server_cur][0]
            list_paths.append(list_paths[-1] + [h])
            server_cur = h
        return list_paths

    def fix_point_backlog_constraints(self, first_server: int, group_flows: List[int], max_burst: float, file):
        """
        Additional constraints for a fixed point when cutting at edges. This uses a stronger result from[Bou19].

        **Warning:** Should be carefully handled to ensure no that this is not used improperly.
        Example of use-case: the ring network, without shaping or greedy-shapers. Here implemented for the
        iterative method.

        :param first_server: the first server of the flows
        :param group_flows: the set of flows for which to write the constraints
        :param file: file where to write the constraints
        :param max_burst: maximum backlog of the grouped flows
        :return: nothing
        """
        file.write('\n/* Constraints for the fixed-point (group)*/\n')
        e = self.next_foi
        list_paths = self._paths_from_server(first_server, group_flows)
        print(list_paths)
        rho_tot = sum([self.network.flows[i].arrival_curve[0].rho for i in group_flows])
        for path in list_paths:
            if path[-1] == self.network.num_servers - 1:
                path_bis = path + [self.network.num_servers]
            else:
                path_bis = path + [self.network.successors[path[-1]][0]]
            for h2 in path_bis[1:]:
                for h1 in path_bis:
                    if h1 < h2:
                        file.write('0 ')
                        for i in group_flows:
                            file.write('+ f{0}s{1}t{2}e{4} - f{0}s{1}t{3}e{4} '.format(i, first_server, h2, h1, e))
                        file.write('<= {0} + {1}t{2}e{4} - {1}t{3}e{4};\n'.format(max_burst, rho_tot, h2, h1, e))
                        for i in group_flows:
                            file.write('+ f{0}s{1}t{2}e{4} - f{0}s{1}t{3}e{4} >= {5}t{2}e{4} - {5}t{3}e{4};\n '
                                       .format(i, first_server, h2, h1, e, self.network.flows[i].arrival_curve[0].rho))

    def write_constraints(self, file):
        """
        Writes all the constraints of the flow of interest.

        :param file: file to write the constraints
        :return: nothing
        """
        file.write('\n/* flow {} */\n'.format(self.list_flows[self.foi]))
        self.time_constraints(file)
        self.arrival_constraints(file)
        self.service_constraints(file)
        self.monotony_constraints(file)
        self.causality_constraints(file)
        self.greedy_shaping_constraints(file)
        self.arrival_shaping_constraints(file)
        self.fix_point_constraints(file)
