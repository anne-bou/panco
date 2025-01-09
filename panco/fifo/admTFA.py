#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/panco

from __future__ import annotations

__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@huawei.com"
__copyright__ = "Copyright (C) 2024, Huawei Technologies France"
__license__ = "BSD-3"

import numpy as np
import subprocess as sp
from typing import List, Tuple
from collections import defaultdict
import re


from panco.descriptor.network import Network
from panco.descriptor.server import Server
from panco.descriptor.flow import Flow
from panco.descriptor.curves import TokenBucket, RateLatency
from panco.descriptor.server import NO, PFR, IR

from panco.lpSolvePath import LPSOLVEPATH


"""
This file is the implementation of the results of the paper:
[Bou] Anne Bouillard, Admission Shaping With Network Calculus. IEEE Netw. Lett. 6(2): 115-118 (2024)
It generalizes the TFA analysis of a FIFO network by 
- allowing shaping flows together when they enter a system
- computing a fixed point analysis when the arrival curves are a minimum of token-bucket curves
- introducing the interleaved and per flow shaping in the resolution of the fixed point, as define in 
[Thomas, Le Boudec, Mifdaoui, RTSS 19].    
"""


def extract_keys(dic, i, j):
    """
    Function that extracts from a dictionary whose entries are triples, the second component, given that the first and
    third are respectively i and j
    :param dic: the dictionary
    :param i: the value of the first component
    :param j: the value of the third component
    :return: the sorted list of the second components

    >>> dictionary = defaultdict()
    >>> dictionary[(0, 1, 2)] = 1
    >>> dictionary[(1, 2, 3)] = 4
    >>> dictionary[(0, 3, 2)] =  0
    >>> extract_keys(dictionary, 0, 2)
    [(1, 1), (3, 0)]
    """
    subset_keys = [b for (a, b, c) in dic.keys() if a == i and c == j]
    li = [(k, dic[(i, k, j)]) for k in subset_keys]
    return sorted(li)


def expand_shaping_fnc(network: Network) -> List[Tuple[int, int, int, Tuple[List[int], int, List[TokenBucket]]]]:
    """
    Expands the shaping at arrival to all servers: if flows [0, 1, 2, 3, 4] are shaped by some arrival curve, and flows
    [0, 1, 2] follow the same initial path until some other server, flows [0, 1, 2] will be shaped at this server, by \
    an arrival curve that can be computed. Here we use the same arrival curves, the burst parameters will be computed \
    in the linear program.
    :param network: the network for the shaping
    :return: the list of shaping and propagation. Each tuple is (number of the shaping, number of the previous shaping \
    (from which it has been built), server of the previous shaping, (server, list of flows, arrival curve))

    >>> flow1 = Flow([TokenBucket(1, 1)], [0, 1])
    >>> flow2 = Flow([TokenBucket(1, 1)], [0, 2])
    >>> flow3 = Flow([TokenBucket(1, 1)], [1, 2])
    >>> flows = [flow1, flow1, flow1, flow2, flow2, flow3, flow3, flow3]
    >>> arr_shaping = (0, [0, 1, 2, 3, 4], [TokenBucket(0, 5)])
    >>> arr_shaping2 = (1, [5, 6, 7], [TokenBucket(0, 3)])
    >>> server = Server([RateLatency(10, 1)], [], 0)
    >>> network1 = Network([server, server, server], flows, [arr_shaping, arr_shaping2])
    >>> expand_shaping_fnc(network1)
    [(2, 0, 0, (1, [0, 1, 2], [0 + 5t])), (3, 0, 0, (2, [3, 4], [0 + 5t])), (4, 1, 1, (2, [5, 6, 7], [0 + 3t]))]
    """
    list_aux1 = list(enumerate(network.arrival_shaping))
    list_aux = [(a, 0, 0, b) for (a, b) in list_aux1]
    list_shaping = []
    c = len(network.arrival_shaping) - 1
    while list_aux:
        (a, b, h, (j, list_flows, shaping)) = list_aux[0]
        for h in network.successors[j]:
            new_flows = list(set(list_flows).intersection(set(network.edges[j, h])))
            if new_flows:
                c += 1
                list_aux += [(c, a, j, (h, new_flows, shaping))]
                list_shaping += [(c, a, j, (h, new_flows, shaping))]
        list_aux = list_aux[1:]
    return list_shaping


class AdmTfa:
    """
    The class TfaLP computes delay bounds using the TFA++ method, using a linear program. For feed-forward networks,
    does the same as [Mifdaoui, Leydier RTSS17], except that it makes no assumption on the maximum service rate.
    For cyclic networks, it computes the same as [Thomas, Le Boudec, Mifdaoui, RTSS 19], without this same assumption,
    and directly computes the fix-point, without iterating.
    However, this takes into account only the first token-bucket of the arrival curves
    :param network: the network to analyse.
    :param filename: name of the file to write the linear program
    expand_shaping is buit from the expand_shaping_function
    """
    def __init__(self, network: Network, filename="tfa_new.lp"):
        self.network = network
        self.filename = filename
        self._ff_equiv = None
        self._solve = None
        self.expand_shaping = expand_shaping_fnc(self.network)

    def regulator_on_path(self, i: int) -> List[int]:
        """
        For flow i, computes the list of the last regulator on its path: this gives the parameters of the regulator \
        taken for the analysis. For example, if the path is [0, 1, 2, 3, 4], and regulators are
        (NO, IR, NO, PFR, NO), \
        then the result is [0, 1, 1, 3, 2, 5], and the regulator for server 1 is according to the arrival curve at \
        server 1 amd the regulator for server 3 is according to the arrival curve at server 2.
        :param i: The flow number.
        :return: the list of the last regulators.

        >>> flow =  Flow([TokenBucket(1, 1)], [0, 1, 2, 3, 4])
        >>> server_no = Server([RateLatency(10, 0)], [], NO)
        >>> server_fpr = Server([RateLatency(10, 0)], [], PFR)
        >>> server_ir = Server([RateLatency(10, 0)], [], IR)
        >>> servers = [server_no, server_ir, server_no, server_fpr, server_no]
        >>> network1 = Network(servers, [flow])
        >>> AdmTfa(network1).regulator_on_path(0)
        [0, 1, 1, 3, 2, 5]
        """

        path = self.network.flows[i].path + [self.network.num_servers]
        list_ac = [path[0]]
        last_regulator = path[0]
        path_length = len(path)
        h = 0
        while h < path_length - 1:
            if self.network.servers[path[h]].regulator == NO:
                list_ac += [path[h + 1]]
            elif self.network.servers[path[h]].regulator == PFR:
                list_ac += [last_regulator]
                last_regulator = path[h+1]
            else:
                list_ac += [path[h]]
                last_regulator = path[h+1]
            h += 1
        return list_ac

    def tfa_variables(self, file):
        """
        Writes the constraints of the fix-point variables (the burst transmitted for each flow to the next server)
        :param file: the file in which the constraints are written
        :return: None
        """

        file.write('\n /* sigma variables*/\n')
        for i in range(self.network.num_flows):
            for (l, j) in enumerate(self.network.path[i]):
                if j == self.network.path[i][0]:
                    for k, tb in enumerate(self.network.flows[i].arrival_curve):
                        file.write('x{0}l{1}s{2} = {3};\n'.format(i, k, j, tb.sigma))
                else:
                    for k, tb in enumerate(self.network.flows[i].arrival_curve):
                        file.write('x{0}l{1}s{2} = x{0}l{1}s{3} + {4} d{3};\n'.format(i, k, j,
                                                                                      self.network.path[i][l - 1],
                                                                                      tb.rho))

    def tfa_shaping_constraints(self, file):
        """
        Writes the shaping constraints, for the extended shaping.
        :param file: file where to write the linear constraints
        :return: None
        """
        file.write('\n/* shaping constraints */\n')
        for (a, (j, list_flows, shaping)) in enumerate(self.network.arrival_shaping):
            for (u, tb) in enumerate(shaping):
                file.write('y{0}l{2}s{1} = {3};\n'.format(a, j, u, tb.sigma))
                for i in list_flows:
                    file.write('+ f{0}s{1}u{1}'.format(i, j))
                file.write('<= y{0}l{2}s{1} + {3}u{1};\n'.format(a, j, u, tb.rho))
        file.write('\n/* propagation of the shaping constraints */\n')
        for (a, b, h, (j, list_flows, shaping)) in self.expand_shaping:
            for (u, tb) in enumerate(shaping):
                file.write('y{0}l{2}s{1} = y{5}l{2}s{3} + {4} d{3};\n'.format(a, j, u, h, tb.rho, b))
                for i in list_flows:
                    file.write('+ f{0}s{1}u{1}'.format(i, j))
                file.write('<= y{0}l{2}s{1} + {3}u{1};\n'.format(a, j, u, tb.rho))

    def tfa_constraints_server(self, file):
        """
        Writes the TFA constraints for each server
        :param file: the file where the constraints are written
        :return: None
        """
        for j in range(self.network.num_servers):
            file.write('\n /* server {0}*/\n'.format(j))
            for i in self.network.flows_in_server[j]:
                for k, tb in enumerate(self.network.flows[i].arrival_curve):
                    # print(type(tb))
                    file.write('f{0}s{1}u{1} <= x{0}l{2}s{1} + {3} u{1};\n'.format(i, j, k, tb.rho))
            for h in self.network.predecessors[j]:
                for tb in self.network.servers[h].max_service_curve:
                    for i in self.network.edges[(h, j)]:
                        file.write('+ f{0}s{1}u{1}'.format(i, j))
                    file.write('<= {0} + {1} u{2};\n'.format(tb.sigma, tb.rho, j))
            file.write("0")
            for i in self.network.flows_in_server[j]:
                file.write('+ f{0}s{1}u{1}'.format(i, j))
            file.write('= a{0}u{0};\n'.format(j))
            for rl in self.network.servers[j].service_curve:
                file.write('b{0}t{0} >= {1} t{0} - {2};\n'.format(j, rl.rate, rl.rate * rl.latency))
            file.write('b{0}t{0} >= 0;\n'.format(j))
            file.write('b{0}t{0} = a{0}u{0};\n'.format(j))
            file.write('d{0} = t{0} - u{0};\n'.format(j))
            file.write('d{0} >= 0;\n'.format(j))

    def regulator_constraints(self, file):
        """
        Constraints for the regulators. Instead of propagating the bursts, when there is a regulator, takes the \
        output process of a previous server according to the regulation.
        :param file: file where to write the constraints
        :return: None
        """
        file.write("\n/* Constraint from the regulator */ \n")
        for i in range(self.network.num_flows):
            reg = self.regulator_on_path(i)
            for (k, j) in enumerate(self.network.flows[i].path):
                if k > 0 and self.network.servers[self.network.flows[i].path[k - 1]].regulator > 0:
                    file.write("x{0}l0s{1} = x{0}l0s{2};\n".format(i, j, reg[k]))

    @property
    def solve(self) -> Tuple[np.ndarray, defaultdict, defaultdict]:
        """
        Computes the delay bounds of all the servers, and the burst parameters for all flows and shaping at all servers.
        :return: the list of the delays of the servers, and of the bursts
        """
        if self._solve is not None:
            return self._solve
        file = open(self.filename, 'w')
        file.write('max:')
        for i in range(self.network.num_servers):
            file.write('+ d{} '.format(i))
        file.write(';\n')
        self.tfa_constraints_server(file)
        self.tfa_variables(file)
        self.tfa_shaping_constraints(file)
        self.regulator_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S2", self.filename],
                   stdout=sp.PIPE, encoding='utf-8').stdout
        tab_values = s.split('\n')[4:-1]
        values = [[token for token in line.split(' ') if not token == ""] for line in tab_values]
        if not values:
            return self.network.num_servers * [np.inf], defaultdict(), defaultdict()
        tab_delays = np.zeros(self.network.num_servers)
        dict_bursts = defaultdict()
        dict_bursts_shaping = defaultdict()
        for [s1, s2] in values:
            if s1[0] == 'd':
                tab_delays[int(float(s1[1:]))] = float(s2)
            if s1[0] == 'x':
                s = [int(a) for a in re.split("[xsl]", s1)[1:]]
                dict_bursts[(s[0], s[1], s[2])] = float(s2)
            if s1[0] == 'y':
                s = [int(a) for a in re.split("[ysl]", s1)[1:]]
                dict_bursts_shaping[(s[0], s[1], s[2])] = float(s2)
        self._solve = tab_delays, dict_bursts, dict_bursts_shaping
        return self._solve

    @property
    def solve_optim(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the delay bounds of all the servers, and the level at which the
        worst-case delay is computed for each \
        server. (Used for optimization)
        :return: the list of the delays of the servers and level at which the worst-case delay happens.
        """
        # print('solve', [self.network.servers[j].regulator for j in range(self.network.num_servers)])
        file = open(self.filename, 'w')
        file.write('max:')
        for i in range(self.network.num_servers):
            file.write('+ d{} '.format(i))
        file.write(';\n')
        self.tfa_constraints_server(file)
        self.tfa_variables(file)
        self.tfa_shaping_constraints(file)
        self.regulator_constraints(file)
        file.close()
        s = sp.run(["wsl", "lp_solve", "-S2", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        tab_values = s.split('\n')[4:-1]
        values = [[token for token in line.split(' ') if not token == ""] for line in tab_values]
        if not values:
            return self.network.num_servers * [np.inf], self.network.num_servers * [np.inf]
        tab_delays = np.zeros(self.network.num_servers)
        tab_x = np.zeros(self.network.num_servers)
        for [s1, s2] in values:
            if s1[0] == 'd':
                tab_delays[int(float(s1[1:]))] = float(s2)
            if s1[0] == 'a':
                s = [int(a) for a in re.split("[au]", s1)[1:]]
                tab_x[s[0]] = float(s2)
        return tab_delays, tab_x

    @property
    def delay_servers(self) -> List[float]:
        """
        Computes the list of delays for all servers
        :return: the list of delays for all servers
        """
        return list(self.solve[0])

    @property
    def dict_bursts(self):
        """
        Computes the dictionary of the bursts parameters for all flows at each server
        :return: the dictionary of the bursts parameters for all flows at each server
        """
        return self.solve[1]

    @property
    def dict_bursts_shaping(self):
        """
        Computes the dictionary of the bursts parameters for all shaping.
        :return: the dictionary of the bursts parameters for all shaping
        """
        return self.solve[2]

    def delay(self, foi: int) -> float:
        """
        Returns the delay bound for flow foi
        :param foi: the flow of interest
        :return: the delay bound of foi
        """
        tab_delays = self.delay_servers
        return sum([tab_delays[j] for j in self.network.path[foi]])

    @property
    def all_delays(self) -> List[float]:
        """
        Returns the delay bounds for all the flows
        :return: the list of delay bounds
        """
        tab_delays = self.delay_servers
        return [sum([tab_delays[i] for i in self.network.path[j]]) for j in range(self.network.num_flows)]

    @property
    def idx_equiv(self):
        """
        This is used to compute the equivalent network (made of elementary 1-server networks). Every path is then \
        split into elementary path, and flows are renumbered. In the equivalent network, the flows are then (i, j) \
        for all flows i and all servers j on the path of i. This computes the dictionary flow re-numbering in the \
        equivalent network.
        :return: dictionary of flow numbering of the equivalent network.
        """
        dict_flows = defaultdict()
        k = 0
        for i in range(self.network.num_flows):
            for j in self.network.flows[i].path:
                dict_flows[(i, j)] = k
                k += 1
        return dict_flows

    @property
    def ff_equiv(self) -> Network:
        """
        The equivalent network decomposed into elementary servers
        :return: a network
        """
        if self._ff_equiv is not None:
            return self._ff_equiv
        flows = []
        for i in range(self.network.num_flows):
            for j in self.network.flows[i].path:
                li = extract_keys(self.dict_bursts, i, j)
                new_arrival_curve = [TokenBucket(s, self.network.flows[i].arrival_curve[k].rho)
                                     for (k, s) in li]
                flows += [Flow(new_arrival_curve, [j])]
        idx = self.idx_equiv
        shaping_constraints = self.network.arrival_shaping + [a[3] for a in self.expand_shaping]
        new_arrival_shaping = []
        for (k, (j, fl, sh)) in enumerate(shaping_constraints):
            li = extract_keys(self.dict_bursts_shaping, k, j)
            new_flows = [idx[i, j] for i in fl]
            new_shaping = [TokenBucket(s, sh[l].rho) for (l, s) in li]
            new_arrival_shaping += [(j, new_flows, new_shaping)]
        for (j, h) in self.network.edges.keys():
            if self.network.servers[j].max_service_curve:
                new_flows = [idx[i, h] for i in self.network.edges[(j, h)]]
                new_arrival_shaping += [(h, new_flows, self.network.servers[j].max_service_curve)]
        self._ff_equiv = Network(self.network.servers, flows, new_arrival_shaping)

        return self._ff_equiv

    def single_server_backlog(self, server: int) -> float:
        ff = self.ff_equiv
        file = open(self.filename, 'w')
        file.write('\n /* LP for server {0}*/\n'.format(server))
        file.write('max: a - b;\n')
        for i in ff.flows_in_server[server]:
            for k, tb in enumerate(ff.flows[i].arrival_curve):
                file.write('f{0} <= {1} + {2} u;\n'.format(i, tb.sigma, tb.rho))
        for (j, list_flows, shaping) in ff.arrival_shaping:
            if j == server:
                file.write("0")
                for tb in shaping:
                    for i in list_flows:
                        file.write('+ f{0}'.format(i))
                    file.write('<= {0} + {1} u;\n'.format(tb.sigma, tb.rho))
        file.write("0")
        for i in ff.flows_in_server[server]:
            file.write(" + f{0}".format(i))
        file.write(" = a;\n")
        for rl in ff.servers[server].service_curve:
            file.write('b >= {1} u - {2};\n'.format(server, rl.rate, rl.rate * rl.latency))
        file.write('b >= 0;\n')
        file.close()
        s = sp.run(LPSOLVEPATH + ["-s5", "-S1", self.filename],
                   stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])

    def single_server_backlogged_period(self, server: int) -> float:
        ff = self.ff_equiv
        file = open(self.filename, 'w')
        file.write('\n /* LP for server {0}*/\n'.format(server))
        file.write('max: t;\n')
        for i in ff.flows_in_server[server]:
            for k, tb in enumerate(ff.flows[i].arrival_curve):
                file.write('f{0} <= {1} + {2} t;\n'.format(i, tb.sigma, tb.rho))
        for (j, list_flows, shaping) in ff.arrival_shaping:
            if j == server:
                file.write("0")
                for tb in shaping:
                    for i in list_flows:
                        file.write('+ f{0}'.format(i))
                    file.write('<= {0} + {1} t;\n'.format(tb.sigma, tb.rho))
        file.write("0")
        for i in ff.flows_in_server[server]:
            file.write(" + f{0}".format(i))
        file.write(" = a;\n")
        for rl in ff.servers[server].service_curve:
            file.write('b >= {1} t - {2};\n'.format(server, rl.rate, rl.rate * rl.latency))
        file.write('b >= 0;\n')
        file.write('b <= a;\n')
        file.close()
        s = sp.run(LPSOLVEPATH + ["-s5", "-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])

    @property
    def all_backlogs(self):
        return [self.single_server_backlog(j) for j in range(self.network.num_servers)]

    @property
    def all_backlogged_periods(self):
        return [self.single_server_backlogged_period(j) for j in range(self.network.num_servers)]

    @property
    def extract_regulators(self):
        """
        Extracts the shaping regulators after ff_equiv hav been computed.
        - for IR, this is the arrival curve at the server
        - for PFR, this is the output of the previous regulator on the path (of the arrival curve)
        :return:the list of the regulators
        """
        ff = self.ff_equiv
        tab_regulators = []
        idx = self.idx_equiv
        for i in range(self.network.num_flows):
            reg_path = self.regulator_on_path(i)
            for (k, j) in enumerate(self.network.flows[i].path):
                if self.network.servers[j].regulator > 0:
                    tab_regulators += [(i, j, self.network.servers[j].regulator,
                                        ff.flows[idx[(i, reg_path[k+1])]].arrival_curve[0])]
        return tab_regulators


if __name__ == '__main__':
    import doctest
    doctest.testmod()
