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


from typing import List

from panco.descriptor.curves import TokenBucket
from panco.descriptor.server import Server
from panco.descriptor.network import Network
from panco.tsn.tsnFlow import TsnFlow, aggregate_tsn_flows, tsn_num_classes
from panco.tsn.tsnSwitch import TsnSwitch

"""
This module decomposes a network made o TSN switches into one FIFO network per traffic
class. 
"""
class TsnNetwork:

    """
    Class that describes a TSN network and transforms in into a list of FIFO networks (one per CBS/BE class)

    :param self.num_servers: number of servers in the network. It corresponds to the number of TSN switches
    :param self.num_classes: number of classes in the network
    :param self.flows: list of TSN flows circulating in the network
    :param self.servers: list of TSN switches (max_packet can be intialized to 0 at the construction step)
    """
    def __init__(self, flows: List[TsnFlow], servers: List[TsnSwitch], be_flows=False, be_max_length=0):
        self.num_servers = len(servers)
        self.flows = flows
        self.servers = servers
        self.be_max_length = be_max_length
        if be_flows:
            self.num_classes = tsn_num_classes(flows)
        else:
            self.num_classes = tsn_num_classes(flows) + 1
        self.be_max_length = be_max_length

    def to_fifo_networks(self) -> List[Network]:
        """
        Function that transforms the TSN network into a list of FIFO networks, one per CBS/BE class. It first updates \\
        the maximum packet length of the flows for each servers in oder to compute residual service curves.

        :return: the list of Network for the analysis of each class.
        """
        list_flows, packet_length = aggregate_tsn_flows(self.flows, self.num_servers,
                                                        self.num_classes, self.be_max_length)
        list_servers = [[] for __ in range(self.num_classes)]
        for j in range(self.num_servers):
            server = self.servers[j]
            self.servers[j] = TsnSwitch(server.bandwidth, server.period, server.tas_intervals,
                                        server.guard_intervals, server.idleSlopes, packet_length[j])
        for j in range(self.num_servers):
            residual_servers = self.servers[j].residual_cbs + [self.servers[j].best_effort_ssc]
            max_shapers = self.servers[j].shaping_cbs + [TokenBucket(0,  self.servers[j].bandwidth)]
            for c in range(self.num_classes):
                list_servers[c] += [Server([residual_servers[c]], [max_shapers[c]])]
        return [Network(list_servers[c], list_flows[c]) for c in range(self.num_classes)]
