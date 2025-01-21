#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Augmentation of exposure-response prevention with transcranial direct current stimulation for contamination-related OCD：a randomized clinical trial
# Time    : 2021-10-10
# Author  : Wenjun Jia
# File    : microstate.py


import numpy as np

class MicrostateParam:
    def __init__(self, sequence, n_microstate):
        self.sequence = sequence
        self.n_microstate = n_microstate
        self.n_sequence = len(sequence)

    def calculate_duration(self, window=None, step=None):
        def duration(sequence):
            res = []
            label = sequence[0]
            count = 1
            res_temp = {}
            for i in range(self.n_microstate):
                res_temp[i] = []
            for j in range(1, len(sequence)):
                if label == sequence[j]:
                    count += 1
                else:
                    label = sequence[j]
                    res_temp[sequence[j - 1]].append(count)
                    count = 1
            for i in range(self.n_microstate):
                res.append(np.sum(res_temp[i]) / len(res_temp[i]) if len(res_temp[i]) > 0 else 1)
            return res

        if window and step:
            n_step = 0
            res = [[] for _ in range(self.n_microstate)]
            while (window+step*n_step) <= self.n_sequence:
                temp = duration(self.sequence[step*n_step:window+step*n_step])
                for index, element in enumerate(temp):
                    res[index].append(element)
                n_step += 1
            return res
        else:
            return duration(self.sequence)

    def calculate_frequency(self, window=None, step=None):
        if window and step:
            n_step = 0
            res = [[] for _ in range(self.n_microstate)]
            while (window+step*n_step) <= self.n_sequence:
                seq = self.sequence[step * n_step:window + step * n_step]
                temp = np.zeros(self.n_microstate)
                i = 0
                j = 0
                while i != len(seq) -1 and j != len(seq)-1:
                    for j in range(i, len(seq)):
                        if seq[i] != seq[j]:
                            temp[seq[i]] += 1
                            i = j
                            break
                temp[seq[j]] += 1
                for i in range(self.n_microstate):
                    res[i].append(temp[i])
                n_step += 1
        else:
            res = []
            res_temp = {}
            for i in range(self.n_microstate):
                res_temp[i] = []
            n_block = int(self.n_sequence / window)
            for i in range(n_block):
                label = self.sequence[i * window]
                temp = {}
                for j in range(i * window + 1, (i + 1) * window):
                    if label != self.sequence[j]:
                        if label in temp:
                            temp[label] += 1
                        else:
                            temp[label] = 1
                        label = self.sequence[j]
                for key, value in temp.items():
                    res_temp[key].append(value)
            for i in range(self.n_microstate):
                res.append(np.mean(res_temp[i]))
        return res

    def calculate_coverage(self, window=None, step=None):
        if window and step:
            n_step = 0
            res = [[] for _ in range(self.n_microstate)]
            while (window+step*n_step) <= self.n_sequence:
                for i in range(self.n_microstate):
                    res[i].append(np.argwhere(np.asarray(self.sequence[step*n_step:window+step*n_step]) == i).shape[0] / window)
                n_step += 1
        else:
            res = []
            for i in range(self.n_microstate):
                res.append(np.argwhere(np.asarray(self.sequence) == i).shape[0] / self.n_sequence)


        return res