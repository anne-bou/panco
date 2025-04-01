##!/usr/bin/env python3
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
from panco.descriptor.curves import TokenBucket, RateLatency


NO = 0
PFR = 1
IR = 2


class Server:
    """
    The Server class encodes the description of a server, characterized by:

        - a (minimal) service curve :math:`\\beta(t) = \\max_i(R_i(t-T_i)_+)`
        - a greedy shaping curve :math:`\\gamma(t) = \\min_i (\\sigma_i + \\rho_i t)`


    :param service_curve: the service curve, given by a maximum of the rate-latency functions :math:`\\beta`.
    :type service_curve: List[RateLatency]
    :param max_service_curve: the maximum curve, given by a minimum of the token-bucket functions :math:`\\gamma`
    :type max_service_curve: List[TokenBucket]
    :param regulator: is the server followed by a regulator? (0=no, 1=per-flow regulator, 2=interleaved regulator).
    Default=0
    :type regulator: int


    >>> Server([RateLatency(10, 1), RateLatency(20, 2)], [TokenBucket(5, 20), TokenBucket(0, 30)])
    <Server: β(t) = max [10(t - 1)_+, 20(t - 2)_+]
             \u03C3(t) = min [5 + 20t, 0 + 30t]>
    <BLANKLINE>
    """
    def __init__(self, service_curve: List[RateLatency], max_service_curve: List[TokenBucket], regulator=0):
        self.service_curve = service_curve
        self.max_service_curve = max_service_curve
        self.regulator = regulator

    def __str__(self) -> str:
        return "β(t) = max %s\n         \u03C3(t) = min %s, reg %s" % (self.service_curve, self.max_service_curve,
                                                                       self.regulator)

    def __repr__(self) -> str:
        return "<Server: %s>\n" % self.__str__()

    def __eq__(self, other: Server):
        return self.max_service_curve == other.max_service_curve and self.service_curve == other.service_curve \
               and self.regulator == other.regulator

    @property
    def is_constant_rate(self) -> bool:
        if not (self.regulator == 0):
            return False
        if not len(self.service_curve) == 1 or not len(self.max_service_curve) == 1:
            return False
        if not self.service_curve[0].latency == 0 and not self.max_service_curve[0].sigma == 0:
            return False
        if self.service_curve[0].rate == self.max_service_curve[0].rho:
            return True
        return False


class ConstantRateServer(Server):
    def __init__(self, rate: float):
        service_curve = [RateLatency(rate, 0)]
        max_service_curve = [TokenBucket(0, rate)]
        super().__init__(service_curve, max_service_curve, NO)
