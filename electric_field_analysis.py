#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Augmentation of exposure-response prevention with transcranial direct current stimulation for contamination-related OCD：a randomized clinical trial
# Time    : 2021-10-10
# Author  : Wenjun Jia
# File    : microstate.py

import os
import numpy as np
import simnibs
import nibabel as nib

class EfAnalysis:
    def __init__(self, subject, gm_surf, m2m_folder):
        self.gm_surf = simnibs.read_msh(gm_surf)
        self.m2m_folder = m2m_folder
        self.subject = subject

    def get_roi_ef_fsaverage_by_area(self, field_name='E_normal', rois=None):
        atlas = simnibs.subject_atlas('HCP_MMP1', self.m2m_folder)
        flag = True
        for sub_rois in rois:
            if sub_rois[3] == 'L':
                roi_name = str(sub_rois[4])
                lh_roi_name = 'lh.' + roi_name
                rh_roi_name = 'rh.' + roi_name
                lh_roi = atlas[lh_roi_name]
                rh_roi = atlas[rh_roi_name]
                if flag:
                    roi = np.logical_or(lh_roi, rh_roi)
                    flag = False
                else:
                    item = np.logical_or(lh_roi, rh_roi)
                    roi = np.logical_or(item, roi)
        self.gm_surf.add_node_field(roi, 'ROI')
        node_areas = self.gm_surf.nodes_areas()
        # ef = np.average(self.gm_surf.field[field_name][roi], weights=node_areas[roi])
        ef = np.median(self.gm_surf.field[field_name][roi])
        return ef

    def get_roi_ef_fsaverage(self, field_name='E_normal', hemisphere=None, roi_name=None):
        atlas = simnibs.subject_atlas('HCP_MMP1', self.m2m_folder)
        if hemisphere:
            roi_name = hemisphere + '.' + roi_name
            roi = atlas[roi_name]
        else:
            lh_roi_name = 'lh.' + roi_name
            rh_roi_name = 'rh.' + roi_name
            lh_roi = atlas[lh_roi_name]
            rh_roi = atlas[rh_roi_name]
            roi = np.logical_or(lh_roi, rh_roi)

        self.gm_surf.add_node_field(roi, 'ROI')
        node_areas = self.gm_surf.nodes_areas()
        ef = np.average(self.gm_surf.field[field_name][roi], weights=node_areas[roi])
        return ef

    def get_roi_ef_individual(self, lh_hcp_annotation_path, rh_hcp_annotation_path, hcp_label_path, hemisphere, field_name):
        lh_vertices = nib.freesurfer.read_annot(lh_hcp_annotation_path)[0]
        rh_vertices = nib.freesurfer.read_annot(rh_hcp_annotation_path)[0]
        lh_vertices = np.zeros(len(lh_vertices), dtype=bool)
        rh_vertices = np.zeros(len(rh_vertices), dtype=bool)
        index, _ = fMRIRoi.get_atlas_position(hcp_label_path)
        if hemisphere == 'lh':
            lh_vertices[index] = True
        elif hemisphere == 'rh':
            rh_vertices[index] = True
        roi = np.concatenate((lh_vertices, rh_vertices))
        self.gm_surf.add_node_field(roi, 'ROI')
        node_areas = self.gm_surf.nodes_areas()
        ef = np.average(self.gm_surf.field[field_name][roi], weights=node_areas[roi])
        return ef
