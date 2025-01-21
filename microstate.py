#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Augmentation of exposure-response prevention with transcranial direct current stimulation for contamination-related OCD：a randomized clinical trial
# Time    : 2021-10-10
# Author  : Wenjun Jia
# File    : microstate.py

import mat73
import numpy as np
import mne
import multiprocessing
from multiprocessing import Pool
from scipy.io import savemat, loadmat


class Microstate:
    def __init__(self, data):
        self.data = Microstate.substract_mean(data)
        self.n_t = self.data.shape[0]
        self.n_ch = self.data.shape[1]
        self.gfp = np.std(self.data, axis=1)
        self.gfp_peaks = Microstate.locmax(self.gfp)
        self.gfp_values = self.gfp[self.gfp_peaks]
        self.n_gfp = self.gfp_peaks.shape[0]
        self.sum_gfp2 = np.sum(self.gfp_values**2)

    @staticmethod
    def locmax(x):
        dx = np.diff(x)
        zc = np.diff(np.sign(dx))
        m = 1 + np.where(zc == -2)[0]
        return m

    @staticmethod
    def substract_mean(x):
        return x - x.mean(axis=1, keepdims=True)

    @staticmethod
    def assign_labels_kmeans(data, maps, n_ch, gfp, gfp_peaks=None):
        c = np.dot(data, maps.T)
        if isinstance(gfp_peaks, np.ndarray):
            c /= (n_ch * np.outer(gfp[gfp_peaks], np.std(maps, axis=1)))
        else:
            c /= (n_ch * np.outer(gfp, np.std(maps, axis=1)))
        l = np.argmax(c ** 2, axis=1)
        return l, c

    def fit_back(self, maps, threshold=None):
        c = np.dot(self.data, maps.T) / (self.n_ch * np.outer(self.gfp, maps.std(axis=1)))
        c = abs(c)
        c_max = np.max(c, axis=1)
        c_max_index = np.argmax(c, axis=1)
        if threshold:
            c_threshold_index = np.where(c_max > threshold)[0]
            l = c_max_index[c_threshold_index]
        else:
            l = c_max_index
        return l

    def fit_back_peaks(self, maps, threshold=None):
        c = np.dot(self.data[self.gfp_peaks], maps.T) / (
                    self.n_ch * np.outer(self.gfp[self.gfp_peaks], maps.std(axis=1)))
        c = abs(c)
        c_max = np.max(c, axis=1)
        c_max_index = np.argmax(c, axis=1)
        res = []
        if threshold:
            c_threshold_index = np.where(c_max > threshold)[0]
            l = c_max_index[c_threshold_index]
            gfp_peaks = self.gfp_peaks[c_threshold_index]
        else:
            l = c_max_index
            gfp_peaks = self.gfp_peaks
        for i in range(0, len(gfp_peaks) - 1):
            med_point = (gfp_peaks[i] + gfp_peaks[i + 1]) // 2
            res += [l[i] for j in range(med_point - gfp_peaks[i])] + [l[i + 1] for j in
                                                                      range(gfp_peaks[i + 1] - med_point)]
        return np.asarray(res)

    def gev(self, maps):
        n_maps = len(maps)
        c = np.dot(self.data[self.gfp_peaks], maps.T)
        c /= (self.n_ch * np.outer(self.gfp[self.gfp_peaks], np.std(maps, axis=1)))
        l = np.argmax(c ** 2, axis=1)
        gev = np.zeros(n_maps)
        for k in range(n_maps):
            r = l == k
            gev[k] = np.sum(self.gfp_values[r] ** 2 * c[r, k] ** 2) / self.sum_gfp2
        return gev, np.sum(gev)

    def wcss(self, maps):
        c = np.dot(self.data[self.gfp_peaks], maps.T) / (self.n_ch * np.outer(self.gfp[self.gfp_peaks], maps.std(axis=1)))
        c = abs(c)
        c_max = np.max(c, axis=1)
        c_max_index = np.argmax(c, axis=1)
        l = c_max_index
        gfp_peaks = self.gfp_peaks


    def kl_criterion(self, wcss):
        kl_values = []
        for k in range(1, len(wcss)-1):
            numerator = np.abs(wcss[k+1] - 2 * wcss[k] + wcss[k-1])
            denominator = wcss[k]
            kl_value = numerator / denominator
            kl_values.append(kl_value)
        return np.array(kl_values)



    def kmeans(self, n_maps, maxerr=1e-6, maxiter=1000):
        np.random.seed()
        rndi = np.random.permutation(self.n_gfp)[:n_maps]
        data_gfp = self.data[self.gfp_peaks, :]
        sum_v2 = np.sum(data_gfp ** 2)
        maps = data_gfp[rndi, :]
        maps /= np.sqrt(np.sum(maps ** 2, axis=1, keepdims=True))
        n_iter = 0
        var0 = 1.0
        var1 = 0.0
        while ((np.abs((var0 - var1) / var0) > maxerr) & (n_iter < maxiter)):
            l_peaks, c = Microstate.assign_labels_kmeans(data_gfp, maps, self.n_ch, self.gfp, self.gfp_peaks)
            for k in range(n_maps):
                vt = data_gfp[l_peaks == k, :]
                sk = np.dot(vt.T, vt)
                evals, evecs = np.linalg.eig(sk)
                v = np.real(evecs[:, np.argmax(np.abs(evals))])
                maps[k, :] = v / np.sqrt(np.sum(v ** 2))
            var1 = var0
            var0 = sum_v2 - np.sum(np.sum(maps[l_peaks, :] * data_gfp, axis=1) ** 2)
            var0 /= (self.n_gfp * (self.n_ch - 1))
            n_iter += 1
        l, _ = Microstate.assign_labels_kmeans(self.data, maps, self.n_ch, self.gfp)
        var = np.sum(self.data ** 2) - np.sum(np.sum(maps[l, :] * self.data, axis=1) ** 2)
        var /= (self.n_t * (self.n_ch - 1))
        cv = var * (self.n_ch - 1) ** 2 / (self.n_ch - n_maps - 1.) ** 2
        return maps, l, cv

    def wrap_kmeans(self, para):
        return self.kmeans(para[0])

    def kmeans_repetition(self, n_repetition, n_maps, n_pool=11):
        l_list = []
        cv_list = []
        maps_list = []
        pool = Pool(n_pool)
        multi_res = []
        for i in range(n_repetition):
            multi_res.append(pool.apply_async(self.wrap_kmeans, ([n_maps],)))
        pool.close()
        pool.join()
        for i in range(n_repetition):
            temp = multi_res[i].get()
            maps_list.append(temp[0])
            l_list.append(temp[1])
            cv_list.append(temp[2])
        k_opt = np.argmin(cv_list)
        return maps_list[k_opt], cv_list[k_opt]

    def microstate(self, max_maps, n_repetition, n_pool=11, is_single=None):
        maps_list = []
        cv_list = []
        if is_single is not None:
            temp = self.kmeans_repetition(n_repetition, is_single, n_pool)
            maps_list.append(temp[0].tolist())
            cv_list.append(temp[1])
        else:
            for n_maps in range(1, max_maps + 1):
                print('n_maps:', n_maps)
                temp = self.kmeans_repetition(n_repetition, n_maps, n_pool)
                maps_list.append(temp[0].tolist())
                cv_list.append(temp[1])
        return maps_list, cv_list


