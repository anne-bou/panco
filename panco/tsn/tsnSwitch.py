#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/panco


__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@huawei.com"
__copyright__ = "Copyright (C) 2024, Huawei Technologies France"
__license__ = "BSD-3"

from __future__ import annotations

from typing import List, Tuple
import numpy as np

from PyLPbounds.Descriptor.curves import TokenBucket, RateLatency

"""
Implementation of a single switch. THe aim is to compute the service curves for the CBS (credit-based scheduler) 
classes. The CBS scheduler is described by an Idle Slope. We also assume that Time-triggered traffic is also 
present, an it is described by a list of gate-opening times. We also assume the 
presence of guard bands.  
"""


def _shift_intervals(period: float, intervals: List[Tuple[float, float]], start: float) -> List[float]:
    """
    Computes a list of intervals when starting at start: shift of the sequence

    :param period: period of the process
    :param intervals: interval list during one period
    :param start: starting date of the interval
    :return: the shifted list of intervals
    """
    if start > period:
        start = np.mod(start, period)
    list_begin = []
    list_end = []
    for (a, b) in intervals:
        if a < start < b:
            raise Exception('start should be an element in intervals')
        if a >= start:
            list_begin += [(a - start, b - start)]
        if b <= start:
            list_end += [(a + period - start, b + period - start)]
    return list_begin + list_end


def _intervals_to_curves(period: float, list_intervals: List[Tuple[float, float]]) -> Tuple[TokenBucket, RateLatency]:
    if not list_intervals == []:
        return TokenBucket(0, 0), RateLatency(0, 0)
    b_max = 0
    lat_max = 0
    # print(list_intervals)
    rho = sum([t2 - t1 for (t1, t2) in list_intervals]) / period
    for (a, b) in list_intervals:
        length_interval = 0
        intervals = _shift_intervals(period, list_intervals, a)
        for (c, d) in intervals:
            length_interval += d - c
            b_max = max(b_max, length_interval - rho * d)
        length_interval = 0
        intervals = _shift_intervals(period, list_intervals, b)
        for (c, d) in intervals:
            lat_max = max(lat_max, c - length_interval / rho)
            length_interval += d - c
    return TokenBucket(b_max, rho), RateLatency(rho, lat_max)


def _contract_time(period: float, guard_intervals: List[Tuple[float, float]],
                   frozen_intervals: List[Tuple[float, float]]) -> Tuple[float, List[Tuple[float, float]]]:
    frozen_time = 0
    contracted_guard = []
    i = 0
    (a_g, b_g) = guard_intervals[0]
    for (a_f, b_f) in frozen_intervals:
        while a_g < a_f:
            if b_g > b_f:
                raise Exception('guard and frozen intervals must be disjoint')
            contracted_guard += [(a_g - frozen_time, b_g - frozen_time)]
            i += 1
            if i < len(guard_intervals):
                (a_g, b_g) = guard_intervals[i]
            else:
                (a_g, b_g) = (np.inf, np.inf)
        frozen_time += b_f - a_f
    return period - frozen_time, contracted_guard


class TsnSwitch:
    """
    Class for a TSN Switch

    :param bandwidth: bandwidth of the switch :math:`R`
    :type bandwidth: float
    :param period: Period of the GCL :math:`P`
    :type period: float
    :param tas_intervals: List of intervals reserved for the TAS (TT) traffic :math:`\\{[s_i, s_i + d_i]\\}`
    :type tas_intervals: List[Tuple[float, float]]
    :param guards_intervals:  List of intervals reserved for the guard bands :math:`\\{[s_i - g_i, s_i]\\}`
    :type guards_intervals: List[Tuple[float, float]]
    :param self.idleSlopes: the IdleSlope parameters of the CBS classes :math:`I_i`
    :type self.idleSlopes: List[float]
    :param max_length: maximum length of a packet of CBS of BE class: for the best effort it is the last element
    :type max_length: List[float]
    :param self.sendSlopes: the SendSlope parameters of the CBS classes. By convention, :math:`S_i = I_i - R`
    :type self.sendSlopes: List[float]
    :param self.num_cbs: number of CBS classes
    :type self.num_cbs: int

    >>> period = 16
    >>> bandwidth = 10
    >>> tas_intervals = [(0, 2), (6, 7), (10, 13)]
    >>> guard_intervals = [(-1, 0), (4.5, 6), (8.5, 10)]
    >>> idleSlopes = [2, 3]
    >>> max_length = [1, 3, 2]
    >>> tsn = TsnSwitch(bandwidth, period, tas_intervals, guard_intervals, idleSlopes, max_length)
    >>> tsn.is_stable
    True
    >>> tsn.tas_load
    0.375
    >>> tsn.tas_curves
    (2.0 + 0.375t, 0.375(t - 5.333333333333334)_+)
    >>> tsn.non_frozen_time_curves
    (0.625(t - 3.2)_+, 2.0 + 0.625t)
    >>> tsn.sendSlopes
    [-8, -7]
    >>> tsn.length_bar
    [3, 2]
    >>> tsn.guard_curves
    (1.2 + 0.4t, 0.4(t - 3.0)_+)
    >>> tsn.min_credit
    [-0.8, -2.1]
    >>> tsn.max_credit
    [5.0, 11.100000000000001]
    >>> tsn.residual_cbs
    [1.25(t - 7.2)_+, 1.875(t - 9.120000000000001)_+]
    >>> tsn.shaping_cbs
    [9.8 + 1.25t, 19.200000000000003 + 1.875t]
    >>> tsn.best_effort_ssc
    3.125(t - 15.68)_+
    """
    def __init__(self, bandwidth: float, period: float, tas_intervals: List[Tuple[float, float]],
                 guards_intervals: List[Tuple[float, float]], idleslopes: List[float], max_length: List[float]):
        self.bandwidth = bandwidth
        self.period = period
        self.tas_intervals = tas_intervals
        self.guard_intervals = guards_intervals
        self.idleSlopes = idleslopes
        self.max_length = max_length
        self.sendSlopes = [ids - self.bandwidth for ids in self.idleSlopes]
        self.num_cbs = len(idleslopes)

    @property
    def tas_load(self) -> float:
        """
        Computes the average amount of time reserved for the TAS/TT traffic: :math:`\\rho_f = \\frac{\\sum_i d_i}{R}`

        :return: :math:`\\rho_f`
        :rtype: float
        """
        sum_tas_periods = sum([t2-t1 for (t1, t2) in self.tas_intervals])
        return sum_tas_periods / self.period

    @property
    def tas_curves(self) -> Tuple[TokenBucket, RateLatency]:
        """
        Computes the period curves for the TAS frozen-credit periods: Finds the optimal pseudo-linear period curves \\
        such that :math:`\\rho_f(t-s-\\tau_f)_+ \\leq F(s, t) \\leq \\sigma_f + \\rho_f (t-s)`.

        :return: :math:`(t\\mapsto \\sigma_f + \\rho_f t, t\\mapsto \\rho_f(t-\\tau_f)_+)
        :rtype: Tuple[TokenBucket, RateLatency]
        """
        return _intervals_to_curves(self.period, self.tas_intervals)

    @property
    def non_frozen_time_curves(self) -> Tuple[RateLatency, TokenBucket]:
        """
        Computes the period curves for the non-TAS frozen-credit periods: Finds the optimal pseudo-linear period \\
        curves such that :math:`\\rho_{nf}(t-s-\\tau_{nf})_+ \\leq F(s, t) \\leq \\sigma_{nf} + \\rho_{nf} (t-s)`.

        :return: :math:`(t\\mapsto \\sigma_{nf} + \\rho_{nf} t, t\\mapsto \\rho_{nf}(t-\\tau_{nf})_+)
        :rtype: Tuple[TokenBucket, RateLatency]
        """
        tas_max, tas_min = self.tas_curves
        nf_min = RateLatency(1 - self.tas_load, tas_max.sigma / (1 - self.tas_load))
        nf_max = TokenBucket(tas_min.latency * self.tas_load, 1 - self.tas_load)
        return nf_min, nf_max

    @property
    def guard_curves(self) -> Tuple[TokenBucket, RateLatency]:
        if not self.guard_intervals:
            return TokenBucket(0, 0), RateLatency(0, 0)
        period_guard, interval_guard = _contract_time(self.period, self.guard_intervals, self.tas_intervals)
        return _intervals_to_curves(period_guard, interval_guard)

    @property
    def is_stable(self) -> bool:
        """
        Checks if the switch is stable (finite credits for the CBS classes), that is if \\
        the credits are bounded: :math:`\\sum_{i=1}^p I_i < R(1-\\rho_g)`

        :return: True if the switch is stable, False otherwise
        """
        guard_load = self.guard_curves[0].rho
        return sum(self.idleSlopes) < self.bandwidth * (1 - guard_load)

    @property
    def length_bar(self) -> List[float]:
        """
        Computes the maximum length of a packet with strictly lower priority than the CBS class of interest \\
        (including the BE traffic): :math:`\\bar{L}_i = \\max (L_j,~ j > i)`

        :return: :math:`\\bar{L}_i`
        :rtype: List[float]
        """
        return [max(self.max_length[i+1:]) for i in range(self.num_cbs)]

    @property
    def min_credit(self) -> List[float]:
        """
        Computes the minimum credit for the CBS classes: \\
        :math:`c_i^{\\min} = \\frac{L_iS_i}{R}`

        :return: :math:`c^{\\min}_i`
        :rtype: List[float]
        """
        return [self.sendSlopes[i] * self.max_length[i] / self.bandwidth for i in range(self.num_cbs)]

    @property
    def max_credit(self) -> List[float]:
        """
        Computes the maximum credit for the CBS classes: \\
        :math:`c_i^{\\max}=\\frac{I_i}{R(1-\\rho_g)-\\sum_{j<i}I_i}(R\\sigma_g + \\bar{L}_i - \\sum_{j<i} {c_j^{min}})`

        :return: :math:`c^{\\max}_i`
        :rtype: List[float]
        """
        if not self.is_stable:
            raise Exception('the switch is not stable, credit is not bounded')
        guard_curve = self.guard_curves[0]
        c_min = self.min_credit
        sum_c_min = [sum([c_min[j] for j in range(i)]) for i in range(self.num_cbs)]
        sum_idle = [sum([self.idleSlopes[j] for j in range(i)]) for i in range(self.num_cbs)]
        return [self.idleSlopes[i] / (self.bandwidth * (1 - guard_curve.rho) - sum_idle[i]) *
                (self.bandwidth * guard_curve.sigma + self.length_bar[i] - sum_c_min[i])
                for i in range(self.num_cbs)]

    @property
    def residual_cbs(self) -> List[RateLatency]:
        """
        Computes the residual service curve for the CBS classes: \\
        :math:`\\beta_i(t) = I_i\\rho_{nf} (t - \\tau_{nf} - \\frac{c^{\\max}_i}{I_i \\rho_{nf}})_+`

        :return: :math:`(\\beta_i)_i`
        :rtype: List[RateLatency]
        """
        non_frozen_curve = self.non_frozen_time_curves[0]
        return [RateLatency(self.idleSlopes[i] * non_frozen_curve.rate, non_frozen_curve.latency +
                            self.max_credit[i] / (self.idleSlopes[i] * non_frozen_curve.rate))
                for i in range(self.num_cbs)]

    @property
    def shaping_cbs(self) -> List[TokenBucket]:
        """
        Computes the shaping curves for the CBS classes: \\
        :math:`\\beta^{sh}_i(t) =  \\rho_{nf}I_it + \\sigma_{nf} I_i +  c_i^{\\max} - c_i^{\\min}`

        :return: :math:`(\\beta^{sh}_i)_i`
        :rtype: List[TokenBucket]
        """
        non_frozen_curve = self.non_frozen_time_curves[1]
        return [TokenBucket(self.max_credit[i] - self.min_credit[i] + self.idleSlopes[i] * non_frozen_curve.sigma,
                            self.idleSlopes[i] * non_frozen_curve.rho) for i in range(self.num_cbs)]

    @property
    def best_effort_ssc(self) -> RateLatency:
        """
        Computes a strict residual service curve for the BE traffic: if shaping cured for the CBS are \\
        :math:`\\beta_i(t) = b_i + r_i t`, then :math:`\\beta_{BE}(t) = R_{BE}(t-T_{BE})_+` with\\
        :math:`R_{BE} = R \\rho_{nf} - \\sum_{i} r_i` and \\
        :math:`T_{BE}= \\frac{R  \\rho_{nf}  \\tau_{nf} + \\sum_i b_i}{R_{BE}}`

        :return: :math:`\\beta_{BE}`
        :rtype: RateLatency
        """
        non_frozen_curve = self.non_frozen_time_curves[0]
        shaping_cbs = self.shaping_cbs
        b = sum([shaping_cbs[i].sigma for i in range(self.num_cbs)])
        r = sum([shaping_cbs[i].rho for i in range(self.num_cbs)])
        rate = self.bandwidth * non_frozen_curve.rate - r
        latency = (self.bandwidth * non_frozen_curve.rate * non_frozen_curve.latency + b) / rate
        return RateLatency(rate, latency)
