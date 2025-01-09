#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/Panco


from __future__ import annotations

__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@huawei.com"
__copyright__ = "Copyright (C) 2024, Huawei Technologies France"
__license__ = "BSD-3"

from typing import List

from panco.descriptor.curves import TokenBucket, sum_ac_list, RateLatency
from panco.descriptor.server import Server
from panco.descriptor.network import Network
from panco.staticpriorities.spFlow import SpFlow, sp_num_classes, aggregate_sp_flows
from panco.staticpriorities.spServer import SpServer
from panco.fifo.admTFA import AdmTfa


class SpNetwork:
    """
    Class that describes a static priority network and transforms in into a list of FIFO networks (one per CBS/BE class)

    :param self.num_servers: number of servers in the network. It corresponds to the number of TSN switches
    :param self.num_classes: number of classes in the network
    :param self.flows: list of priority flows circulating in the network
    :param self.servers: (list of servers)
    """
    def __init__(self, servers: List[Server], flows: List[SpFlow]):
        self.num_servers = len(servers)
        self.num_classes = sp_num_classes(flows)
        self.flows = flows
        self.servers = servers
        self.per_class_networks = []
        self.per_class_residual_servers = []
        self.top_class_cross_traffic = [[TokenBucket(0, 0)] for _ in range(self.num_servers)]
        self.per_class_flows, self.max_length = aggregate_sp_flows(flows, self.num_servers)
        self.sp_servers = [SpServer(servers[j].service_curve, servers[j].max_service_curve, self.max_length[j])
                           for j in range(self.num_servers)]

    def per_class_network(self, k: int, is_strict: bool) -> Network:
        """
        Function that computes an equivalent network for priority class k. It uses the  cross traffic of the higher
        priority classes (self.top_class_cross_traffic[j]) for each server j, and computes the residual service for all
        the servers.
        @param k: the class number
        @param is_strict: if strict service curves need to be computed (for example if the scheduling policy inside the
        class need  a strict service curve (DRR or  CBS)). If it is a FIFO network inside the class, then this is not
        needed.
        @return: the equivalent network for class k.
        """
        list_servers_class = []
        for j in range(self.num_servers):
            list_servers_class += [self.sp_servers[j].residual(self.top_class_cross_traffic[j], k, is_strict)]
        return Network(list_servers_class, self.per_class_flows[k])

    def equiv_network(self, is_strict):
        """
        Computes the equivalent network for each priority class. It starts fom the highest priority, computes the
        equivalent for th next class, and continues iteratively.
        @param is_strict: f strict service curves need to be computed (for example if the scheduling policy inside the
        class need  a strict service curve (DRR or CBS)). If it is a FIFO network inside the class, then this is not
        needed.
        @return: the list of networks or each class (they can be used for a FIFO analysis for example).

        >>> flow1 = SpFlow([TokenBucket(1, 1)], [0], 1, 0)
        >>> flow2 = SpFlow([TokenBucket(1, 1)], [1], 1, 1)
        >>> flow3 = SpFlow([TokenBucket(2, 2)], [0, 1], 1, 2)
        >>> server1 = Server([RateLatency(5, 0)], [TokenBucket(0, 10)])
        >>> server2 = Server([RateLatency(5, 0)], [TokenBucket(0, 10)])
        >>> sp_network = SpNetwork([server1, server2], [flow1, flow2, flow3])
        >>> sp_network.equiv_network(True)
        [<Network:
        Flows:
              0: α(t) = min [1 + 1t]; π = [0]
        Servers:
              0: β(t) = max [5(t - 0.2)_+]
                 σ(t) = min [0 + 10t]
              1: β(t) = max [5(t - 0.2)_+]
                 σ(t) = min [0 + 10t]>, <Network:
        Flows:
              0: α(t) = min [1 + 1t]; π = [1]
        Servers:
              0: β(t) = max [4(t - 0.5)_+]
                 σ(t) = min [0 + 10t]
              1: β(t) = max [5(t - 0.2)_+]
                 σ(t) = min [0 + 10t]>, <Network:
        Flows:
              0: α(t) = min [2 + 2t]; π = [0, 1]
        Servers:
              0: β(t) = max [4(t - 0.5)_+]
                 σ(t) = min [0 + 10t]
              1: β(t) = max [4(t - 0.5)_+]
                 σ(t) = min [0 + 10t]>]
        """
        e_net = []
        for k in range(self.num_classes):
            net_elem = self.per_class_network(k, is_strict)
            e_net += [net_elem]
            fifo = AdmTfa(e_net[-1])
            equiv = fifo.ff_equiv
            for j in range(equiv.num_servers):
                if equiv.flows_in_server[j]:
                    tb = sum_ac_list([equiv.flows[i].arrival_curve for i in equiv.flows_in_server[j]])
                    self.top_class_cross_traffic[j] = sum_ac_list([self.top_class_cross_traffic[j], tb])
                    if not self.top_class_cross_traffic[j]:
                        self.top_class_cross_traffic[j] += [TokenBucket(0, 0)]
        return e_net


if __name__ == '__main__':
    import doctest
    doctest.testmod()
