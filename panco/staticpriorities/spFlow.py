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

from typing import List, Tuple
import numpy as np

from panco.descriptor.curves import TokenBucket
from panco.descriptor.flow import Flow

"""
Description of a flow with static priority scheduling, described by  
- its path
- its arrival curve 
- its maximum packet length
- its class : the lower the class number, the higher the priority: class 0 is the highest priority, and the best effort 
class (lower priority), the class with the highest number.     
"""


class SpFlow:
    """
    Class describing a Static priority flow

    :param arrival_curve: arrival curve of the flow
    :param path: path (sequence of TSN switches) of the flow
    :param max_length: maximum packet length of the flow
    :param sp_class: priority class of the flow


    """
    def __init__(self, arrival_curve: List[TokenBucket], path: List[int], max_length: float, sp_class: int):

        self.arrival_curve = arrival_curve
        self.path = path
        self.max_length = max_length
        self.sp_class = sp_class


def sp_num_classes(list_flows: List[SpFlow]) -> int:
    """
    From a list of priority  flows computes the number of  classes,

    :param list_flows: list of priority flows
    :return: the number of classes

    >>> flow1 = SpFlow([TokenBucket(1, 2)], [0, 1], 1024, 0)
    >>> flow2 = SpFlow([TokenBucket(2, 1)], [1, 2], 1024, 1)
    >>> flow3 = SpFlow([TokenBucket(1, 3)], [0, 1, 2], 1024, 0)
    >>> sp_num_classes([flow1, flow2, flow3])
    2
    """
    return max([f.sp_class for f in list_flows]) + 1


def aggregate_sp_flows(list_flows: List[SpFlow], num_servers: int) -> \
        Tuple[List[List[Flow]], np.ndarray]:
    """
    From a list of priority flows and number of servers, computes the list of flows per priority
    classe, and the maximum packet size for each class and each priority link.

    :param list_flows: list of priority flows
    :param num_servers: total number of servers
    :return: the list of flows per class, and the maximum packet length

    >>> flow1 = SpFlow([TokenBucket(1, 2)], [0, 1], 1024, 0)
    >>> flow2 = SpFlow([TokenBucket(2, 1)], [1, 2], 508, 1)
    >>> flow3 = SpFlow([TokenBucket(1, 3)], [0, 1, 2], 1024, 0)
    >>> aggregate_sp_flows([flow1, flow2, flow3], 3)
    ([[<Flow: α(t) = min [1 + 2t]; π = [0, 1]>
    , <Flow: α(t) = min [1 + 3t]; π = [0, 1, 2]>
    ], [<Flow: α(t) = min [2 + 1t]; π = [1, 2]>
    ]], array([[1024.,    0.],
           [1024.,  508.],
           [1024.,  508.]]))
    """
    num_classes = sp_num_classes(list_flows)
    list_flows_per_class = [[] for __ in range(num_classes)]
    max_packet_length = np.zeros((num_servers, num_classes))
    for f in list_flows:
        list_flows_per_class[f.sp_class] += [Flow(f.arrival_curve, f.path)]
        for p in f.path:
            max_packet_length[p, f.sp_class] = max(f.max_length, max_packet_length[p, f.sp_class])
    return list_flows_per_class, max_packet_length


if __name__ == '__main__':
    import doctest
    doctest.testmod()
