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

from panco.fifoCplex.admTFA import AdmTfaCpx
from panco.fifoCplex.elpConstraints import ELPConstraintsCpx
from panco.fifoCplex.plpConstraints import PLPConstraintsCpx
from panco.fifo.sfaLP import SfaLP
from panco.fifoCplex.sfaLP import SfaLPCpx
# from panco.fifo.tfaLP import TfaLP
from panco.lpSolvePath import LPSOLVEPATH
from panco.fifo.admTFA import AdmTfa
from panco.fifoCplex.cplex_lex import sol_extraction, lex_cplex, extract_obj, cplex_solve



class TreeLPCpx:
    # Linear analysis for fifo tree networks
    def __init__(self, network, foi, polynomial=True, sfa=False, tfa=False, filename="fifocplex.lp"):
        self.network = network
        self.foi = foi
        # self.constraints = LPConstraints(network, foi)
        if sfa:
            delay_sfa = SfaLPCpx(network).all_delays
        else:
            delay_sfa = None
        if tfa:
            delay_tfa = AdmTfaCpx(network).delay_servers
            # print(max(delay_tfa))
        else:
            delay_tfa = None
        if polynomial:
            self.constraints = PLPConstraintsCpx(network, foi, None, None, delay_sfa, delay_tfa)
        else:
            self.constraints = ELPConstraintsCpx(network, foi)
        self.filename = filename

    def burst_constraints(self, file):
        for i in range(self.network.num_flows):
            file.write('x{0} = {1}\n'.format(i, self.network.flows[i].arrival_curve[0].sigma))

    def delay_objective(self, file):
        if self.network.path[self.foi][-1] == self.network.num_servers - 1:
            file.write('maximize \n obj: t0e0 - t{}e0\n'.format(self.constraints.t_min[self.network.path[self.foi][0]]))
        else:
            file.write('flow do not stop at last server\n')
        file.write('subject to\n')

    @property
    def delay(self):
        file = open(self.filename, 'w')
        self.delay_objective(file)
        self.constraints.time_constraints(file)
        self.constraints.arrival_constraints(file)
        self.constraints.fifo_constraints(file)
        self.constraints.service_constraints(file)
        self.constraints.monotony_constraints(file)
        self.constraints.shaping_constraints(file)
        self.constraints.arrival_shaping_constraints(file, True)
        self.constraints.sfa_delay_constraints(file)
        self.constraints.tfa_delay_constraints(file)
        self.burst_constraints(file)
        file.write('End')
        file.close()
        s = cplex_solve(self.filename)
        filename = 'soly'
        d = extract_obj(filename)
        return d
        # dict_sol = lex_cplex(filename)
        # tab_delays, dict_bursts, dict_bursts_shaping = sol_extraction(dict_sol, self.network)
        # self._solve = tab_delays, dict_bursts, dict_bursts_shaping
        # return(self._solve)
        # s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        # return float(s.split()[-1])

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
        self.constraints.arrival_shaping_constraints(file, True)
        self.constraints.sfa_delay_constraints(file)
        self.constraints.tfa_delay_constraints(file)
        self.burst_constraints(file)
        file.write('End')
        file.close()
        # print(self.filename)
        s = cplex_solve(self.filename)
        filename = 'soly'
        d = extract_obj(filename)
        return d
        # file.close()
        # s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        # return float(s.split()[-1])

    def backlog_set_of_flows(self, set_flows):
        file = open(self.filename, 'w')
        self.constraints.backlog_set_objective(set_flows, file)
        self.constraints.time_constraints(file)
        self.constraints.arrival_constraints(file)
        self.constraints.fifo_constraints(file)
        self.constraints.service_constraints(file)
        self.constraints.monotony_constraints(file)
        self.constraints.shaping_constraints(file)
        self.constraints.arrival_shaping_constraints(file, True)
        self.constraints.sfa_delay_constraints(file)
        self.constraints.tfa_delay_constraints(file)
        self.burst_constraints(file)
        file.write('End')
        file.close()
        s = cplex_solve(self.filename)
        filename = 'soly'
        d = extract_obj(filename)
        return d
        # file.close()
        # s = sp.run(LPSOLVEPATH + ["-S1", self.filename], stdout=sp.PIPE, encoding='utf-8').stdout
        # return float(s.split()[-1])
