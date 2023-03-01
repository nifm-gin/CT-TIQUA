#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 16:29:58 2021

@author: cbrossard
"""

import bids
from nipype.interfaces import fsl

from multiprocessing import Pool
from shutil import copyfile
#from cmtklib.interfaces.fsl import FSLCreateHD
import ants



        
def Apply_ANTS_Atlas(DataPath, Atlas, Template, OutAtlasesPath, OutTemplatePath, n_jobs):

    layout = bids.BIDSLayout(DataPath, validate = False)
    layout_df = layout.to_df()
    sub_df = layout_df[(layout_df.extension == 'nii.gz')]
    subsub_df = sub_df[(sub_df.session == 'J0')]
    
    #ref_path = '/media/cbrossard/ClementBackUp1/Database_CT_Radiomic_Test/Database_CT_Radiomic_Test/BIDS_Radiomic_TBI/derivatives/FLIRT_on_Images_masked/sub-P3_ses-J1_masked_seg.nii.gz'
    #pool = Pool(n_jobs)
    
    
    for index, row in subsub_df.iterrows():
        
        #if (row.subject in {'Patient01', 'Patient02'}) or (row.subject=='Patient05' and row.session == 'D03') or (row.subject=='Patient0' and row.session == 'D01'):
        #if row.subject not in {'P01', 'P02'}:
        #if row.subject not in {'P09'}:
            print(row)
            img_fixed = ants.image_read(row.path)
            img_moving = ants.image_read(Template)
            reg = ants.registration(img_fixed, img_moving)
            outfilename = 'sub-' + row.subject + '_ses-' + row.session + '_TemplateRegistered.nii.gz'
            reg['warpedmovout'].to_file(OutTemplatePath+outfilename)
            
            mytx = reg['fwdtransforms']
            im_to_embarque = ants.image_read(Atlas)
            embarqued_im = ants.apply_transforms(img_fixed, im_to_embarque, transformlist=mytx, interpolator='nearestNeighbor')
            outatlas = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_AtlasRegistered.nii.gz'
            embarqued_im.to_file(outatlas)
            #pool.apply_async(applyxfm.run)

            
            
    #pool.close()
    #pool.join() 


def Apply_ANTS_Atlas_2(DataPath, Atlas_path, Template_path, OutAtlasesPath, OutTemplatePath, n_jobs):

    layout = bids.BIDSLayout(DataPath, validate = False)
    layout_df = layout.to_df()
    sub_df = layout_df[(layout_df.extension == 'nii.gz')]
    subsub_df = sub_df[(sub_df.session == 'J0')]
    
    #ref_path = '/media/cbrossard/ClementBackUp1/Database_CT_Radiomic_Test/Database_CT_Radiomic_Test/BIDS_Radiomic_TBI/derivatives/FLIRT_on_Images_masked/sub-P3_ses-J1_masked_seg.nii.gz'
    #pool = Pool(n_jobs)
    layout_T = bids.BIDSLayout(Template_path, validate = False)
    layout_A = bids.BIDSLayout(Atlas_path, validate = False)
    
    for index, row in subsub_df.iterrows():
        
        #if (row.subject in {'Patient01', 'Patient02'}) or (row.subject=='Patient05' and row.session == 'D03') or (row.subject=='Patient0' and row.session == 'D01'):
        #if row.subject not in {'P01', 'P02'}:
        #if row.subject not in {'P09'}:
            print(row.path)
            if not row.subject == 'P04':
                continue
            Tmpl = layout_T.get(subject = row.subject, session= row.session, return_type='filename', extension='nii.gz')
            img_fixed = ants.image_read(row.path)
            img_moving = ants.image_read(Tmpl[0])
            reg = ants.registration(img_fixed, img_moving)
            outfilename = 'sub-' + row.subject + '_ses-' + row.session + '_TemplateRegistered.nii.gz'
            reg['warpedmovout'].to_file(OutTemplatePath+outfilename)
            
            mytx = reg['fwdtransforms']
            Atl = layout_A.get(subject = row.subject, session= row.session, return_type='filename', extension='nii.gz')
            im_to_embarque = ants.image_read(Atl[0])
            embarqued_im = ants.apply_transforms(img_fixed, im_to_embarque, transformlist=mytx, interpolator='nearestNeighbor')
            outatlas = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_AtlasRegistered.nii.gz'
            embarqued_im.to_file(outatlas)
            #pool.apply_async(applyxfm.run)

            
            
    #pool.close()
    #pool.join() 