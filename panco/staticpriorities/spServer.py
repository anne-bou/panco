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

from panco.descriptor.curves import TokenBucket, RateLatency, residual_general, sum_ac
from panco.descriptor.server import Server

"""
Description of a static priority server by
- a minimum (strict)service curve 
- a maximum service curve (usually a greedy shaper)
- the number of class
- the maximum packet length for each class
"""


class SpServer:
    """
    Class describing a static priority server:
    :param service_curve: the minimal service curve of the server
    :param max_service_curve: the maximum service curve of the server
    :param max_length: the maximum packet length for each static priority class
    """
    def __init__(self, service_curve: List[RateLatency], max_service_curve: List[TokenBucket],
                 max_length: List[float]):
        self.service_curve = service_curve
        self.max_service_curve = max_service_curve
        self.max_length = max_length
        self.num_classes = len(max_length)

    @property
    def length_bar_simple(self) -> List[float]:
        """
        Computes the maximum length of a packet with strictly lower priority. This is used when
        computing a (min,plus) residual service curve.

        :return: the maximum packet length of classes with strictly higher priority, for each class.

        >>> sp_server = SpServer([RateLatency(5, 0)], [TokenBucket(1, 10)], [508, 1024, 508])
        >>> sp_server.length_bar_simple
        [1024, 508, 0]
        """
        return [max(self.max_length[i + 1:]) for i in range(self.num_classes - 1)] + [0]

    @property
    def length_bar_strict(self) -> List[float]:
        """
        Computes the maximum length of a packet with lower or equal priority. This is used when
        a strict residual service curve is omputed

        :return: the maximum packet length of classes with higher or equal priority, for each class.

        >>> sp_server = SpServer([RateLatency(5, 0)], [TokenBucket(1, 10)], [508, 1024, 508])
        >>> sp_server.length_bar_strict
        [1024, 1024, 508]
        """
        return [max(self.max_length[i:]) for i in range(self.num_classes)]

    def residual(self, cross_traffic: List[TokenBucket], i: int, is_strict: bool) -> Server:
        """
        Computes the residual service curve of class i, assuming this is the highest priority,
        in the server (all other priorities have been handled). Two cases are considered: whether
        the residual needs to be strict or not.
        @param cross_traffic: the arrival curve of the cross traffic of class i
        @param i: the class to remove from the service
        @param is_strict: True if the residual service needs to be strict, False otherwise.
        @return: the server with the residual service curve

        >>> sp_server = SpServer([RateLatency(5, 0)], [TokenBucket(1, 10)], [3, 1, 1])
        >>> sp_server.residual([TokenBucket(1, 1)], 0,  True)
        <Server: β(t) = max [4(t - 1.0)_+]
                 σ(t) = min [1 + 10t]>
        <BLANKLINE>
        >>> sp_server = SpServer([RateLatency(5, 0)], [TokenBucket(1, 10)], [3, 1, 1])
        >>> sp_server.residual([TokenBucket(1, 1)], 0,  False)
        <Server: β(t) = max [4(t - 0.5)_+]
                 σ(t) = min [1 + 10t]>
        <BLANKLINE>
        """
        if is_strict:
            ac = sum_ac(cross_traffic, [TokenBucket(self.length_bar_strict[i], 0)])
        else:
            ac = sum_ac(cross_traffic, [TokenBucket(self.length_bar_simple[i], 0)])
        return Server(residual_general(self.service_curve, ac), self.max_service_curve)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
