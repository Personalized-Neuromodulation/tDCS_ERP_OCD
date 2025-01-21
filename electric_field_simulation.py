#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Augmentation of exposure-response prevention with transcranial direct current stimulation for contamination-related OCD：a randomized clinical trial
# Time    : 2021-10-10
# Author  : Wenjun Jia
# File    : microstate.py

import os
import numpy as np
from simnibs import sim_struct, run_simnibs

class EfSimulation:
    def __init__(self, head_msh, pathfem, eeg_cap):
        self.head_msh = head_msh
        self.pathfem = pathfem
        self.eeg_cap = eeg_cap

    def tDCS(self, electrodes, shape='ellipse', dimensions=[10, 10], thickness=4):
        s = sim_struct.SESSION()
        s.fnamehead = self.head_msh
        s.map_to_surf = True
        s.pathfem = self.pathfem
        s.eeg_cap = self.eeg_cap
        s.open_in_gmsh = False
        # s.fields = 'eEjJ'
        tdcslist = s.add_tdcslist()
        tdcslist.currents = electrodes['current']
        for i in range(len(electrodes['current'])):
            ele = tdcslist.add_electrode()
            ele.channelnr = i + 1
            ele.centre = electrodes['name'][i]
            ele.shape = shape
            ele.dimensions = dimensions
            ele.thickness = thickness
        run_simnibs(s)