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
from typing import List

from panco.descriptor.network import Network
from panco.descriptor.flow import Flow
from panco.descriptor.server import Server
from panco.descriptor.curves import TokenBucket, RateLatency
from panco.blind.blindConstraints import BlindConstraints
from panco.lpSolvePath import LPSOLVEPATH

class TreeBlindLP:
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
        self.constraints = BlindConstraints(network, foi)

    def burst_constraints(self, file):
        """
        Write the burst constraints (in a tree, bursts are equal to the sigma parameter of the flows)

        :param file: file where to write the constraints
        :return: nothing
        """
        file.write('\n/* the burst constraints*/\n')
        for i in range(self.network.num_flows):
            file.write('x{0} = {1};\n'.format(i, self.network.flows[i].arrival_curve[0].sigma))

    def delay_objective(self, file):
        """
        Writes the delay objective and some additional constraints for the flow of interest

        :param file: file where constraints are written
        :return: nothing
        """
        if not self.network.path[self.foi][-1] == self.network.num_servers - 1:
            raise Exception('flow do not stop at last server')
        file.write('max: t{0}e0 - ue0;\n'.format(self.network.num_servers))
        file.write('ue0 <= t{0}e0;\n'.format(self.network.num_servers))
        j = self.network.path[self.foi][0]
        file.write('ue0 >= t{0}e0;\n'.format(j))
        file.write('f{0}s{1}t{1}e0 = f{0}s{2}ue0;\n'.format(self.foi, self.network.num_servers, j))
        for tb in self.network.flows[self.foi].arrival_curve:
            file.write('f{0}s{1}ue0 - f{0}s{1}t{1}e0 <= {2} + {3}ue0 - {3}t{1}e0;\n'.
                        format(self.foi, j, tb.sigma, tb.rho))
        file.write('f{0}s{1}ue0 - f{0}s{1}t{1}e0 >= 0;\n'.format(self.foi, j))

    def backlog_objective(self, file):
        """
        Writes the backlog objective of the flow of interest.

        :param file: file where constraints are written
        :return: nothing
        """
        if self.network.path[self.foi][-1] == self.network.num_servers - 1:
            file.write('max: f{0}s{1}t{2}e0 - f{0}s{2}t{2}e0;\n'.format(self.foi, self.network.flows[self.foi].path[0],
                                                                        self.network.num_servers))
        else:
            raise Exception('flow do not stop at last server\n')

    def backlog_set_objective(self, set_flows: List[int], file):
        """
        Writes the backlog objective of a set of flows of interest flow of interest.

        :param set_flows: set of flows of interest
        :param file: the file where constraints are written
        :return: nothing
        """
        file.write('max: ')
        for i in set_flows:
            if self.network.path[i][-1] == self.network.num_servers - 1:
                file.write('+ f{0}s{1}t{2}e0 - f{0}s{2}t{2}e0'.format(i, self.network.flows[i].path[0],
                                                                     self.network.num_servers))
            else:
                raise Exception('flow do not stop at last server\n')
        file.write(';\n')


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
        >>> TreeBlindLP(tree, 0, 'test_tree_blind.lp').delay
        10.06361149
        """
        file = open(self.filename, 'w')
        self.delay_objective(file)
        self.constraints.time_constraints(file)
        self.constraints.arrival_constraints(file)
        self.constraints.causality_constraints(file)
        self.constraints.service_constraints(file)
        self.constraints.greedy_shaping_constraints(file)
        self.constraints.monotony_constraints(file)
        self.constraints.arrival_shaping_constraints(file)
        self.burst_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
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
        >>> TreeBlindLP(tree, 0, 'test_tree_blind_backlog.lp').backlog
        11.03419973
        """
        file = open(self.filename, 'w')
        self.backlog_objective(file)
        self.constraints.time_constraints(file)
        self.constraints.arrival_constraints(file)
        self.constraints.causality_constraints(file)
        self.constraints.service_constraints(file)
        self.constraints.greedy_shaping_constraints(file)
        self.constraints.monotony_constraints(file)
        self.constraints.arrival_shaping_constraints(file)
        self.burst_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        # s = sp.run(["wsl", "lp_solve", "-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])

    def backlog_set_of_flows(self, set_flows: List[int], max_burst: float) -> float:
        """
        Computes the worst-case backlog of a set of flows of interest.

        **Warning:** do not use if you do not master [Bou19], implemented for the iterative method.

        :param set_flows: set of flows of interest
        :type set_flows: List[int]
        :param max_burst: maximum burst during the current iteration
        :type max_burst: float
        :return: the maximum backlog of the set of flows for the next iteration
        :rtype: float
        """
        file = open(self.filename, 'w')
        self.backlog_set_objective(set_flows, file)
        self.constraints.time_constraints(file)
        self.constraints.arrival_constraints(file)
        self.constraints.causality_constraints(file)
        self.constraints.service_constraints(file)
        self.constraints.greedy_shaping_constraints(file)
        self.constraints.monotony_constraints(file)
        self.constraints.fix_point_backlog_constraints(0, [i + 1 for i in set_flows], max_burst, file)
        self.burst_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        # s = sp.run(["wsl", "lp_solve", "-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])


if __name__ == '__main__':
    import doctest
    doctest.testmod()
