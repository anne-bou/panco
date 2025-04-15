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
import subprocess as sp

from panco.descriptor.network import Network
from panco.descriptor.flow import Flow
from panco.descriptor.server import Server
from panco.descriptor.curves import TokenBucket, RateLatency
from panco.lpSolvePath import LPSOLVEPATH


class TreeOnlyBlindLP:
    """
    Class for the analysis for tree networks with blind multiplexing using linear programming techniques.

    :param network: the network to analyze (must be a well-numbered tree)
    :type network: Network
    :param foi: the flow of interest for computing the performance bound
    :type foi: int
    :param filename: the name of the file where the linear program is written
    :type filename: str
    """
    def __init__(self, network: Network, foi: int, filename="tree_blind.lp"):
        """
        Constructor for the class TreeBlindLP
        :param network: the network to analyze (must be a well-numbered tree)
        :param foi: the flow of interest for computing the performance bound
        :param filename: the name of the file where the linear program is written
        """
        self.network = network
        self.foi = foi
        self.filename = filename

    def time_constraints(self, f):
        """
        Writes the time constraints for the network.

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Time Constraints */\n')
        f.write('t{0} <= t{1};\n'.format(self.network.num_servers - 1, self.network.num_servers))
        for j in range(self.network.num_servers - 1):
            h = self.network.successors[j][0]
            f.write('t{0} <= t{1};\n'.format(j, h))

    def arrival_constraints(self, f):
        """
        Writes the arrival constraints for each flow of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* arrival constraints */\n')
        for i in range(self.network.num_flows):
            j = self.network.path[i][0]
            for tb in self.network.flows[i].arrival_curve:
                for (n, h1) in enumerate(self.network.path[i]):
                    if not self.network.successors[self.network.path[i][-1]]:
                        sub_p = self.network.path[i][n + 1:] + [self.network.num_servers]
                    else:
                        sub_p = self.network.path[i][n + 1:] + [self.network.successors[self.network.path[i][-1]][0]]
                    for h2 in sub_p:
                        f.write('f{0}s{1}t{2} - f{0}s{1}t{3} <= {4} + {5} t{2} - {5} t{3};\n'.
                                format(i, j, h2, h1, tb.sigma, tb.rho))

    def arrival_shaping_constraints(self, f):
        """
        Writes the arrival shaping constraints for each flow of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* arrival shaping constraints */\n')
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
                            f.write('+ f{0}s{1}t{2} - f{0}s{1}t{3} '.format(i, j, h2, h1))
                        f.write('<= {0} + {1}t{2} - {1}t{3};\n'.format(tb.sigma, tb.rho, h2, h1))

    def monotony_constraints(self, f):
        """
        Writes the monotony constraints for each arrival process of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Monotony constraints */\n')
        for i in range(self.network.num_flows):
            j = self.network.path[i][0]
            h = self.network.path[i][-1]
            if h == self.network.num_servers - 1:
                s = self.network.num_servers
            else:
                s = self.network.successors[h][0]
            path = self.network.path[i] + [s]
            for (n, h) in enumerate(self.network.path[i]):
                f.write('f{0}s{1}t{2} - f{0}s{1}t{3} <= 0; \n'.format(i, j, h, path[n + 1]))

    def causality_constraints(self, f):
        """
        Writes the causality constraints for each arrival process of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Causality constraints */\n')
        for i in range(self.network.num_flows):
            j = self.network.path[i][0]
            h = self.network.path[i][-1]
            if h == self.network.num_servers - 1:
                s = self.network.num_servers
            else:
                s = self.network.successors[h][0]
            path = self.network.path[i][1:] + [s]
            for h in path:
                f.write('f{0}s{1}t{3} - f{0}s{2}t{3} >= 0; \n'.format(i, j, h, h))

    def service_constraints(self, f):
        """
        Writes the service constraints for each server of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Service constraints */\n')
        for j in range(self.network.num_servers):
            if j == self.network.num_servers - 1:
                h = self.network.num_servers
            else:
                h = self.network.successors[j][0]
            for rl in self.network.servers[j].service_curve:
                for i in self.network.flows_in_server[j]:
                    f.write('f{0}s{1}t{1} - f{0}s{2}t{2} + '.format(i, h, j))
                f.write('{0} >= {1} t{2} - {1} t{3};\n'.format(rl.rate * rl.latency, rl.rate, h, j))
                for i in self.network.flows_in_server[j]:
                    f.write('f{0}s{1}t{1} >= f{0}s{2}t{2};\n '.format(i, h, j))

    def greedy_shaping_constraints(self, f):
        """
        Writes the greedy-shaping constraints for each server of the network

        :param f: file to write the constraints
        :return: nothing
        """
        f.write('\n/* Greedy Shaping constraints */\n')
        for j in range(self.network.num_servers - 1):
            h = self.network.successors[j][0]
            for tk in self.network.servers[j].max_service_curve:
                f.write('0')
                for i in self.network.edges[(j, h)]:  # flows_in_server[j]:
                        f.write('+ f{0}s{1}t{1} - f{0}s{2}t{2}'.format(i, h, j))
                f.write('<= {3} + {0} t{1} - {0} t{2};\n'.format(tk.rho, h, j, tk.sigma))

    def write_constraints(self, file):
        """
        Writes all the constraints of the flow of interest.

        :param file: file to write the constraints
        :return: nothing
        """
        self.time_constraints(file)
        self.arrival_constraints(file)
        self.service_constraints(file)
        self.monotony_constraints(file)
        self.causality_constraints(file)
        self.greedy_shaping_constraints(file)
        self.arrival_shaping_constraints(file)

    def delay_objective(self, file):
        """
        Writes the delay objective and some additional constraints for the flow of interest

        :param file: file where constraints are written
        :return: nothing
        """
        if not self.network.path[self.foi][-1] == self.network.num_servers - 1:
            raise Exception('flow do not stop at last server')
        file.write('max: t{0} - u;\n'.format(self.network.num_servers))
        file.write('u <= t{0};\n'.format(self.network.num_servers))
        j = self.network.path[self.foi][0]
        file.write('u >= t{0};\n'.format(j))
        file.write('f{0}s{1}t{1} = f{0}s{2}u;\n'.format(self.foi, self.network.num_servers, j))
        for tb in self.network.flows[self.foi].arrival_curve:
            file.write('f{0}s{1}u - f{0}s{1}t{1} <= {2} + {3}u - {3}t{1};\n'.
                        format(self.foi, j, tb.sigma, tb.rho))
        file.write('f{0}s{1}u - f{0}s{1}t{1} >= 0;\n'.format(self.foi, j))

    def backlog_objective(self, file):
        """
        Writes the backlog objective of the flow of interest.

        :param file: file where constraints are written
        :return: nothing
        """
        if self.network.path[self.foi][-1] == self.network.num_servers - 1:
            file.write('max: f{0}s{1}t{2} - f{0}s{2}t{2};\n'.format(self.foi, self.network.flows[self.foi].path[0],
                                                                        self.network.num_servers))
        else:
            raise Exception('flow do not stop at last server\n')

    def time_horizon_objective(self, file):
        """
        Writes the time-horizon for the flow of interest

        :param file: file where constraints are written
        :return: nothing
        """
        if not self.network.path[self.foi][-1] == self.network.num_servers - 1:
        #     raise Exception('flow do not stop at last server')
            h = self.network.successors[self.network.path[self.foi][-1]][0]
            file.write('max: t{0} - t{1};\n'.format(h, self.network.path[self.foi][0]))
        else:
            file.write('max: t{0} - t{1};\n'.format(self.network.num_servers, self.network.path[self.foi][0]))
        # file.write('u <= t{0};\n'.format(self.network.num_servers))
        # j = self.network.path[self.foi][0]
        # file.write('u >= t{0};\n'.format(j))
        # file.write('f{0}s{1}t{1} = f{0}s{2}u;\n'.format(self.foi, self.network.num_servers, j))
        # for tb in self.network.flows[self.foi].arrival_curve:
        #    file.write('f{0}s{1}u - f{0}s{1}t{1} <= {2} + {3}u - {3}t{1};\n'.
        #                format(self.foi, j, tb.sigma, tb.rho))
        #file.write('f{0}s{1}u - f{0}s{1}t{1} >= 0;\n'.format(self.foi, j))

    @property
    def delay(self) -> float:
        """
        Computes the worst-case delay of the flow of interest

        **WARNING:** Flows must have a single token-bucket. otherwise, put them in the arrival shaping.

        :return: the delay of the flow of interest
        :rtype: float

        >>> flows = [Flow([TokenBucket(1, 1)], [3, 4]), Flow([TokenBucket(2, 2)], [0, 3]),
        ...          Flow([TokenBucket(3, 3)], [2, 4]), Flow([TokenBucket(4, 4)], [1, 3, 4])]
        >>> servers = [Server([RateLatency(10, 1), RateLatency(5, 0)], []), Server([RateLatency(20, 2)], []),
        ...            Server([RateLatency(30, 3)], []), Server([RateLatency(40, 3)], [TokenBucket(1, 40)]),
        ...            Server([RateLatency(50, 5)], [])]
        >>> arrival_shaping = [(3, [0], [TokenBucket(0, 10)]), (1, [3], [TokenBucket(0, 20)]),
        ...                    (2, [2], [TokenBucket(0, 30)]), (0, [1], [TokenBucket(0, 40)])]
        >>> tree = Network(servers, flows, arrival_shaping)
        >>> TreeOnlyBlindLP(tree, 0, 'test_tree_blind.lp').delay
        10.06361149
        """
        file = open(self.filename, 'w')
        self.delay_objective(file)
        self.time_constraints(file)
        self.arrival_constraints(file)
        self.service_constraints(file)
        self.monotony_constraints(file)
        self.causality_constraints(file)
        self.greedy_shaping_constraints(file)
        self.arrival_shaping_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename],
                   stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])
        # s = sp.run(["wsl", "lp_solve", "-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        # return float(s.split()[-1])

    @property
    def backlog(self) -> float:
        """
        Computes the worst-case backlog of the flow of interest

        :return: the backlog of the flow of interest
        :rtype: float


        >>> flows = [Flow([TokenBucket(1, 1)], [3, 4]), Flow([TokenBucket(2, 2)], [0, 3]),
        ...          Flow([TokenBucket(3, 3)], [2, 4]), Flow([TokenBucket(4, 4)], [1, 3, 4])]
        >>> servers = [Server([RateLatency(10, 1), RateLatency(5, 0)], []), Server([RateLatency(20, 2)], []),
        ...            Server([RateLatency(30, 3)], []), Server([RateLatency(40, 3)], [TokenBucket(1, 40)]),
        ...            Server([RateLatency(50, 5)], [])]
        >>> arrival_shaping = [(3, [0], [TokenBucket(0, 10)]), (1, [3], [TokenBucket(0, 20)]),
        ...                    (2, [2], [TokenBucket(0, 30)]), (0, [1], [TokenBucket(0, 40)])]
        >>> tree = Network(servers, flows, arrival_shaping)
        >>> TreeOnlyBlindLP(tree, 0, 'test_tree_blind_backlog.lp').backlog
        11.03419973
        """
        file = open(self.filename, 'w')
        self.backlog_objective(file)
        self.time_constraints(file)
        self.arrival_constraints(file)
        self.service_constraints(file)
        self.monotony_constraints(file)
        self.causality_constraints(file)
        self.greedy_shaping_constraints(file)
        self.arrival_shaping_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH +["-S1", self.filename],
                   stdout=sp.PIPE, encoding='utf-8').stdout
        # s = sp.run(["wsl", "lp_solve", "-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])

    @property
    def time_horizon(self) -> float:
        """
        Computes the worst-case delay of the flow of interest

        **WARNING:** Flows must have a single token-bucket. otherwise, put them in the arrival shaping.

        :return: the time-horizon of the flow of interest
        :rtype: float

        """
        file = open(self.filename, 'w')
        self.time_horizon_objective(file)

        self.time_constraints(file)

        self.arrival_constraints(file)
        self.service_constraints(file)
        self.monotony_constraints(file)
        self.causality_constraints(file)
        self.greedy_shaping_constraints(file)
        self.arrival_shaping_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename],
                   stdout=sp.PIPE, encoding='utf-8').stdout
        #return float(s.split()[-1])
        # s = sp.run(["wsl", "lp_solve", "-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])
