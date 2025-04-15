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
from typing import List, Tuple

from panco.descriptor.network import Network
from panco.blind.treeBlindLP import TreeBlindLP
from panco.blind.blindConstraints import BlindConstraints
from panco.descriptor.curves import TokenBucket
from panco.descriptor.flow import Flow


class BlindLP:
    """
    Class for the analysis of a network with the linear programming methods.
    The network is decomposed into a forest (self.forest) with well-numbered servers, regarding the decomposition.

    :param network: the network to analyze
    :type network: Network
    :param list_edges: the edges to keep in the decomposition in a forest
    :type list_edges:  List[Tuple[int, int]]
    :param filename: name of the file to write the linear program
    :type filename: str
    """
    def __init__(self, network: Network, list_edges=None, filename="blind.lp"):
        """
        Constructor for the class BlindLP, for the analysis of a network with the linear programming methods.
        The network is decomposed into a forest (self.forest) with well-numbered servers, regarding the decomposition
        """
        self.network = network
        self.filename = filename
        if list_edges is None:
            self.list_edges = list(self.network.edges.keys())
        else:
            self.list_edges = list_edges
        self.forest, self.list_first, z = self.network.decomposition(self.list_edges)
        self.ff_equiv_computed = False
        self._ff_equiv = Network([], [])
        self.stable = True

    def lp_constraint_flow(self, foi: int, file):  # foi flow of the decomposition
        """
        Writes the linear constraints for flow foi of the forest decomposition

        :param foi: flow of interest
        :param file: file where the constraints are written
        :return: None
        """
        net, new_foi, list_flows, list_severs = self.forest.sub_network(foi)
        lp = BlindConstraints(net, new_foi, foi + 1, list_flows)
        lp.write_constraints(file)

    def lp_constraints(self, file):
        """
        Writes the constraints linear program in file

        :param file: the file where the linear program is written
        :return: None
        """
        file.write('max: ')
        for i in range(self.forest.num_flows):
            file.write('+ x{}'.format(i))
        file.write(';\n\n')
        i = 0
        f = 0
        while i < self.forest.num_flows:
            if i == self.list_first[f]:
                file.write('x{0} = {1};\n'.format(i, self.network.flows[f].arrival_curve[0].sigma))
            else:
                self.lp_constraint_flow(i - 1, file)
            i += 1
            if i in self.list_first:
                f += 1

    @property
    def lp_program(self) -> np.ndarray:
        """
        Writes the linear program and solves it to obtain the unknown burst where the flows have been cut

        :return: the list of bursts of flows in the forest
        """
        file = open(self.filename, 'w')
        self.lp_constraints(file)
        file.close()
        s = sp.run(["wsl", "lp_solve", "-S2", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        tab_values = s.split('\n')[4:-1]
        values = [[token for token in line.split(' ') if not token == ""] for line in tab_values]
        tab_bursts = np.inf * np.ones(self.forest.num_flows)
        for [s1, s2] in values:
            if s1[0] == 'x':
                tab_bursts[int(float(s1[1:]))] = float(s2)
        for i in range(self.forest.num_flows):
            if tab_bursts[i] == np.inf:
                self.stable = False
        return tab_bursts

    def update_sigma(self, f: int, sigma: np.ndarray) -> np.ndarray:
        """
        If the network is feed-forward, computes the burst parameters where the flows are cut. Checks if the burst of \\
        the flow has already been computed, otherwise, generate the linear programs needed for the computation.

        :param f: (cut) flow analyzed, its backlog is computed
        :type f: int
        :param sigma: the current array of backlogs for cur flows
        :type sigma: np.ndarray
        :return: the new array of backlogs for cut flows
        :rtype: np.ndarray
        """
        if not sigma[f] == np.inf:
            return sigma
        sub_net, new_f, list_flows, list_servers = self.forest.sub_network(f - 1)
        for j in list_flows:
            if sigma[j] == np.inf:
                sigma = self.update_sigma(j, sigma)
        sigma[f] = TreeBlindLP(sub_net, new_f).backlog
        self.forest.flows[f].arrival_curve[0].sigma = sigma[f]
        return sigma

    def ff_analysis(self) -> Network:
        """
        Analysis performed if the network is feed-forward.

        :return: the updated version of self.forest (updated burst parameters)
        :rtype: Network
        """
        sigma = np.inf * np.ones(self.forest.num_flows)
        for i in range(self.network.num_flows):
            sigma[self.list_first[i]] = self.network.flows[i].arrival_curve[0].sigma
        i = 0
        for f in range(self.forest.num_flows):
            if i < self.network.num_flows and self.list_first[i] == f:
                i += 1
            else:
                self.update_sigma(f, sigma)
        return self.forest

    @property
    def ff_equiv(self) -> Network:
        """
        Construct the equivalent network by solving the fix-point equations. If the network has not been decomposed,
        then returns the original network

        :return: the equivalent network
        """
        if self.ff_equiv_computed:
            return self._ff_equiv
        if self.forest.num_flows == self.network.num_flows:
            self._ff_equiv = self.network
            self.ff_equiv_computed = True
            return self.network
        if self.network.is_feed_forward:
            self._ff_equiv = self.ff_analysis()
            self.ff_equiv_computed = True
            return self._ff_equiv
        else:
            new_sigma = self.lp_program
        for i in range(self.forest.num_flows):
            if i not in self.list_first:
                self.forest.flows[i].arrival_curve[0].sigma = new_sigma[i]
        self._ff_equiv = self.forest
        self.ff_equiv_computed = True
        return self._ff_equiv

    def delay(self, foi: int) -> float:
        """
        Returns the delay bounds for flows foi

        :param foi: the flow of interest
        :return: the delay bound of foi
        """
        ff = self.ff_equiv
        self.ff_equiv_computed = True
        if not self.stable:
            return np.inf
        i = self.list_first[foi]
        residual_latency = 0
        residual_service_rate = np.inf
        while (foi < self.network.num_flows - 1 and i < self.list_first[foi + 1]) or \
              (foi == self.network.num_flows - 1 and i < ff.num_flows):
            tree, foi1, list_flows, list_servers = ff.sub_network(i)
            tree.flows[foi1] = Flow([TokenBucket(0, 0)], tree.path[foi1])
            residual_latency += TreeBlindLP(tree, foi1, "blind_tree.lp").delay
            residual_service_rate = min(residual_service_rate, tree.residual_rate(foi1))
            i += 1
        return residual_latency + self.network.flows[foi].arrival_curve[0].sigma / residual_service_rate

    @property
    def all_delays(self) -> List[float]:
        """
        Returns the delay bounds for all the flows

        :return: the list of delay bounds
        """
        return [self.delay(i) for i in range(self.network.num_flows)]

    def backlog(self, foi: int) -> float:
        """
        Computes the backlog of the flow of interest

        :param foi: the flow of interest
        :type foi: int
        :return: the backlog of foi
        :rtype: float
        """
        ff = self.ff_equiv
        self.ff_equiv_computed = True
        if foi < self.network.num_flows - 1:
            i = self.list_first[foi + 1] - 1
        else:
            i = ff.num_flows
        tree, foi1, list_flows, list_servers = ff.sub_network(i)
        return TreeBlindLP(tree, foi1, "blind_tree.lp").backlog
