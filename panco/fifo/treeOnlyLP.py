#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/panco

from __future__ import annotations

__author__ = "Anne Bouillard"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@ens.fr"
__copyright__ = "Copyright (C) 2026"
__license__ = "BSD-3"

import subprocess as sp

from panco.fifo.plpConstraintsTree import PLPConstraintsTree
from panco.fifo.elpConstraintsTree import ELPConstraintsTree
from panco.fifo.sfaLP import SfaLP
from panco.fifo.admTFA import AdmTfa
from panco.lpSolvePath import LPSOLVEPATH


class TreeOnlyLP:
    # Linear programming analysis for fifo tree networks
    """
    Class for the analysis for tree networks with FIFO scheduling using linear programming techniques.
    This class should be used only for trees, or unfolded feed-forward networks: cutting the network while allowing the
    non-Token-bucket arrival curves could lead to invalid results. Note that the first publications already allow this,
    and this has not been implemented before only for safety reasons: preventing use of this class with cuts.

    :param network: the network to analyze (must be a well-numbered tree)
    :type network: Network
    :param foi: the flow of interest for computing the performance bound
    :type foi: int
    :param filename: the name of the file where the linear program is written
    :type filename: str
    """
    def __init__(self, network, foi, polynomial=True, sfa=True, tfa=True, filename="fifo.lp"):
        self.network = network
        self.foi = foi
        # self.constraints = LPConstraints(network, foi)
        if sfa:
            delay_sfa = SfaLP(network).all_delays
        else:
            delay_sfa = None
        if tfa:
            delay_tfa = AdmTfa(network).delay_servers
        else:
            delay_tfa = None
        if polynomial:
            self.constraints = PLPConstraintsTree(network, foi, delay_sfa, delay_tfa)
        else:
            self.constraints = ELPConstraintsTree(network, foi)
        self.filename = filename

    # def burst_constraints(self, file):
    #     for i in range(self.network.num_flows):
    #         file.write('x{0} = {1};\n'.format(i, self.network.flows[i].arrival_curve[0].sigma))
    #
    # def delay_objective(self, file):
    #     if self.network.path[self.foi][-1] == self.network.num_servers - 1:
    #         file.write('max: t0e0 - t{}e0;\n'.format(self.constraints.t_min[self.network.path[self.foi][0]]))
    #     else:
            #file.write('flow do not stop at last server\n')

    @property
    def delay(self):
        file = open(self.filename, 'w')
        self.constraints.delay_objective(file)
        self.constraints.time_constraints(file)
        self.constraints.arrival_constraints(file)
        self.constraints.fifo_constraints(file)
        self.constraints.service_constraints(file)
        self.constraints.monotony_constraints(file)
        self.constraints.shaping_constraints(file)
        self.constraints.arrival_shaping_constraints(file)
        self.constraints.sfa_delay_constraints(file)
        self.constraints.tfa_delay_constraints(file)

        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])

    @property
    def backlog(self):
        file = open(self.filename, 'w')
        self.constraints.backlog_objective(file)
        self.constraints.time_constraints(file)
        self.constraints.arrival_constraints(file)
        self.constraints.fifo_constraints(file)
        self.constraints.service_constraints(file)
        self.constraints.monotony_constraints(file)
        self.constraints.shaping_constraints(file)
        self.constraints.arrival_shaping_constraints(file)
        self.constraints.sfa_delay_constraints(file)
        self.constraints.tfa_delay_constraints(file)
        file.close()
        s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        return float(s.split()[-1])


