#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of the panco project.
# https://github.com/anne-bou/panco

__author__ = "Anne Bouillard, Ludovic Thomas (Loria)"
__maintainer__ = "Anne Bouillard"
__email__ = "anne.bouillard@huawei.com"
__copyright__ = "Copyright (C) 2022, Huawei Technologies France"
__license__ = "BSD-3"


import platform
import shutil

def get_lpsolve_path():
    if platform.system() == "Linux":
        if shutil.which("lp_solve"):
            return ["lp_solve", "-s5"]
        else:
            return None
    elif platform.system() == "Windows":
        if shutil.which("wsl"):
            return ["wsl", "lp_solve", "-s5"]
        else:
            return None
    else:
        return None
LPSOLVEPATH = get_lpsolve_path()
