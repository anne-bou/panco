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

from __future__ import annotations

from typing import List, Tuple
import numpy as np

from panco.descriptor.curves import TokenBucket
from panco.descriptor.flow import Flow

"""
Description of a TSN flow, described by  
- its maximum packet length
- its class (CBS or BE: if there are p CBS classes from 0 to p-1, then BE is class p
- its path
- its arrival curve 
"""


class TsnFlow:
    """
    Class describing a TSN flow

    :param arrival_curve: arrival curve of the flow
    :param path: path (sequence of TSN switches) of the flow
    :param max_length: maximum packet length of of flow
    :param tsn_class: CBS/BE class of the flow
    """
    def __init__(self, arrival_curve: List[TokenBucket], path: List[int], max_length: float, tsn_class: int):

        self.arrival_curve = arrival_curve
        self.path = path
        self.max_length = max_length
        self.tsn_class = tsn_class


def tsn_num_classes(list_flows: List[TsnFlow]) -> int:
    """
    From a list of TSN flows computes the number of CBS classes (if there is a BE
    class, it is not counted.

    :param list_flows: list of tsn flows
    :return:the number of TSN classes
    """
    return max([f.tsn_class for f in list_flows]) + 1


def aggregate_tsn_flows(list_flows: List[TsnFlow], num_servers: int, num_classes: int, be_max_length: float) -> \
        Tuple[List[List[Flow]], np.ndarray]:
    """
    From a list of TSN flows and number of servers, computes the list of flows per CBS/BE classes, and the maximum \\
    packet size for each class and each tsn link.

    :param list_flows: list o TSN flows
    :param num_servers: total number of servers
    :param be_max_length: maximum packet length of the best_effort traffic
    :return: the list of flows per class, and the maximum packet length
    """
    # num_classes = tsn_num_classes(list_flows)
    list_flows_per_class = [[] for __ in range(num_classes)]
    max_packet_length = np.zeros((num_servers, num_classes))
    for j in range(num_servers):
        max_packet_length[j, num_classes - 1] = be_max_length
    for f in list_flows:
        list_flows_per_class[f.tsn_class] += [Flow(f.arrival_curve, f.path)]
        for p in f.path:
            max_packet_length[p, f.tsn_class] = max(f.max_length, max_packet_length[p, f.tsn_class])
    return list_flows_per_class, max_packet_length
