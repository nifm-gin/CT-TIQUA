#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 19:59:31 2020

@author: cbrossard
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 18:05:46 2020

@author: cbrossard
"""

import bids
from nipype.interfaces import fsl

from multiprocessing import Pool
from shutil import copyfile
#from cmtklib.interfaces.fsl import FSLCreateHD




def Apply_FLIRT(DataPath, ROIPath, OutDataPath, OutROIPath, n_jobs, OtherPath="", OutOtherPath = "", Ref_path=""):

    layout = bids.BIDSLayout(DataPath, validate = False)
    layout_df = layout.to_df()
    sub_df = layout_df[(layout_df.extension == 'nii.gz')]
    
    
    ROIlayout = bids.BIDSLayout(ROIPath, validate = False)
    ROIlayout_df = ROIlayout.to_df()
    ROI_sub_df = ROIlayout_df[(ROIlayout_df.extension == 'nii.gz')]
    
    Otherlayout = bids.BIDSLayout(OtherPath, validate = False)
    Otherlayout_df = Otherlayout.to_df()
    Other_sub_df = Otherlayout_df[(Otherlayout_df.extension == 'nii')]
    
    if Ref_path != "":
        Reflayout = bids.BIDSLayout(Ref_path, validate = False)
        Reflayout_df = Reflayout.to_df()
        Ref_sub_df = Reflayout_df[(Reflayout_df.extension == 'nii.gz')]
    
    paths = sub_df.path
    #ref_path =paths[sub_df.index[0]]  #index du premier fichier.
    #ref_path = paths[sub_df[(sub_df.subject=='Patient01') & (sub_df.session=='D00')].index[0]]
    #ref_path = '/media/cbrossard/ClementBackUp1/Database_CT_Radiomic_Test/Database_CT_Radiomic_Test/BIDS_Radiomic_TBI3/derivatives/Masked_Images/sub-P1_ses-J0_CraneSansIV-masked.nii.gz'
    #ref_path = '/media/cbrossard/ClementBackUp1/Database_CT_Radiomic_Test/Database_CT_Radiomic_Test/BIDS_Radiomic_TBI/derivatives/FLIRT_on_Images_masked/sub-P3_ses-J1_masked_seg.nii.gz'
    pool = Pool(n_jobs)
    
    
    for index, row in sub_df.iterrows():
        #if row.subject =="P26" and row.session=='J3':
            flt = fsl.FLIRT()
            flt.inputs.in_file = row.path
            #ref_path = layout.get(subject = row.subject, session= 'J0', return_type='filename', extension='nii.gz')
            ref_path = Reflayout.get(subject = row.subject, session= row.session, return_type='filename', extension='nii.gz')
            flt.inputs.reference = ref_path[0]
            outfilename = 'sub-' + row.subject + '_ses-' + row.session + '_' + row.suffix + '_registered.nii.gz'
            flt.inputs.out_file = OutDataPath+outfilename
            flt.inputs.out_matrix_file = OutROIPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-matrix.mat'
            flt.inputs.dof = 12
            flt.inputs.bins = 256
            #flt.inputs.cost_func = 'leastsq'
            flt.inputs.cost_func = 'normcorr'
            flt.inputs.interp = 'nearestneighbour'
            #flt.inputs.interp = 'trilinear'
            flt.inputs.searchr_x = [-180, 180]
            flt.inputs.searchr_y = [-180, 180]
            flt.inputs.searchr_z = [-180, 180]
            copyfile(row.path.replace(".nii.gz", ".json"),OutDataPath+outfilename.replace(".nii.gz", ".json"))
            pool.apply_async(flt.run)
    pool.close()
    pool.join() 
            
    pool = Pool(n_jobs)
    for index, row in ROI_sub_df.iterrows():
        #if row.subject =="P25" and row.session=='J3':
        #if (row.subject in {'Patient01', 'Patient02'}) or (row.subject=='Patient05' and row.session == 'D03') or (row.subject=='Patient06' and row.session == 'D01'):
            applyxfm = fsl.ApplyXFM()
            applyxfm.inputs.in_matrix_file = OutROIPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-matrix.mat'
            applyxfm.inputs.in_file = row.path
            applyxfm.inputs.out_file = OutROIPath+'sub-' + row.subject + '_ses-' + row.session + '_' + row.suffix + '_registered.nii.gz'
            #ref_path = ROIlayout.get(subject = row.subject, session= 'J0', return_type='filename', extension='nii.gz')
            ref_path = Reflayout.get(subject = row.subject, session= row.session, return_type='filename', extension='nii.gz')
            applyxfm.inputs.reference = ref_path[0]
            #applyxfm.inputs.reference = row.path
            applyxfm.inputs.apply_xfm = True
            applyxfm.inputs.uses_qform=False
            applyxfm.inputs.no_resample = False
            applyxfm.inputs.interp = 'nearestneighbour'
            pool.apply_async(applyxfm.run)
        
        
    pool.close()
    pool.join()
    if OtherPath != "" and OutOtherPath != "":
        pool = Pool(n_jobs)
        for index, row in Other_sub_df.iterrows():
            #if row.subject =="P25" and row.session=='J3':
            #if (row.subject in {'Patient01', 'Patient02'}) or (row.subject=='Patient05' and row.session == 'D03') or (row.subject=='Patient06' and row.session == 'D01'):
                applyxfm = fsl.ApplyXFM()
                applyxfm.inputs.in_matrix_file = OutROIPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-matrix.mat'
                applyxfm.inputs.in_file = row.path
                applyxfm.inputs.out_file = OutOtherPath+'sub-' + row.subject + '_ses-' + row.session + '_' + row.suffix + '_registered.nii.gz'
                #ref_path = Otherlayout.get(subject = row.subject, session= 'J0', return_type='filename', extension='nii.gz')
                ref_path = Reflayout.get(subject = row.subject, session= row.session, return_type='filename', extension='nii.gz')
                applyxfm.inputs.reference = ref_path[0]
                #applyxfm.inputs.reference = row.path
                applyxfm.inputs.apply_xfm = True
                applyxfm.inputs.uses_qform=False
                applyxfm.inputs.no_resample = False
                applyxfm.inputs.interp = 'nearestneighbour'
                #applyxfm.run()
                pool.apply_async(applyxfm.run)
            
            
        pool.close()
        pool.join()
        
        
        
        
        

        
def Apply_FLIRT_Atlas(DataPath, Atlas, Template, OutAtlasesPath, OutTemplatePath, n_jobs):

    layout = bids.BIDSLayout(DataPath, validate = False)
    layout_df = layout.to_df()
    sub_df = layout_df[(layout_df.extension == 'nii.gz')]
    subsub_df = sub_df[(sub_df.session == 'J0')]
    
    #ref_path = '/media/cbrossard/ClementBackUp1/Database_CT_Radiomic_Test/Database_CT_Radiomic_Test/BIDS_Radiomic_TBI/derivatives/FLIRT_on_Images_masked/sub-P3_ses-J1_masked_seg.nii.gz'
    pool = Pool(n_jobs)
    
    
    for index, row in subsub_df.iterrows():
        #if (row.subject in {'Patient01', 'Patient02'}) or (row.subject=='Patient05' and row.session == 'D03') or (row.subject=='Patient0' and row.session == 'D01'):
        #if row.subject not in {'P01', 'P02'}:
        #if row.subject in {'P09'}:
            flt = fsl.FLIRT()
            print(row)
            flt.inputs.in_file = Template
            flt.inputs.reference = row.path
            outfilename = 'sub-' + row.subject + '_ses-' + row.session + '_TemplateRegistered.nii.gz'
            flt.inputs.out_file = OutTemplatePath+outfilename
            flt.inputs.out_matrix_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-matrix.mat'
            flt.inputs.dof = 7
            flt.inputs.bins = 256
            #flt.inputs.cost_func = 'leastsq'
            flt.inputs.cost_func = 'normcorr'
            flt.inputs.interp = 'nearestneighbour'
            flt.inputs.searchr_x = [-180, 180]
            flt.inputs.searchr_y = [-180, 180]
            flt.inputs.searchr_z = [-180, 180]
            #copyfile(row.path.replace(".nii.gz", ".json"),OutDataPath+outfilename.replace(".nii.gz", ".json"))
            pool.apply_async(flt.run)
            
            
    pool.close()
    pool.join() 
            
    pool = Pool(n_jobs)
    
    for index, row in subsub_df.iterrows():
        #if row.subject not in {'P01', 'P02'}:
        #if row.subject in {'P09'}:
            print(row)
            applyxfm = fsl.ApplyXFM()
            applyxfm.inputs.in_matrix_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-matrix.mat'
            applyxfm.inputs.in_file = Atlas
            applyxfm.inputs.out_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_AtlasRegistered.nii.gz'
            #ref_path = ROIlayout.get(subject = row.subject, session= 'J0', return_type='filename', extension='nii.gz')
            applyxfm.inputs.reference = row.path
            applyxfm.inputs.apply_xfm = True
            applyxfm.inputs.interp = 'nearestneighbour'
            pool.apply_async(applyxfm.run)
            #applyxfm.run()
            
            
    pool.close()
    pool.join() 
    
    
# def Apply_FNIRT_Atlas(DataPath, Atlas, Template, OutAtlasesPath, OutTemplatePath, n_jobs):

#     layout = bids.BIDSLayout(DataPath, validate = False)
#     layout_df = layout.to_df()
#     sub_df = layout_df[(layout_df.extension == 'nii.gz')]
#     subsub_df = sub_df[(sub_df.session == 'J0')]
    
#     #ref_path = '/media/cbrossard/ClementBackUp1/Database_CT_Radiomic_Test/Database_CT_Radiomic_Test/BIDS_Radiomic_TBI/derivatives/FLIRT_on_Images_masked/sub-P3_ses-J1_masked_seg.nii.gz'
#     # pool = Pool(n_jobs)
    
    
#     for index, row in subsub_df.iterrows():
#         #if (row.subject in {'Patient01', 'Patient02'}) or (row.subject=='Patient05' and row.session == 'D03') or (row.subject=='Patient0' and row.session == 'D01'):
#         if row.subject == 'P01':
#             fnt = fsl.FNIRT()
#             fnt.inputs.in_file = Template
#             fnt.inputs.ref_file = row.path
#             fnt.inputs.warped_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-field.nii.gz'
#             fnt.run()
#             #pool.apply_async(fnt.run)
            
            
#     # pool.close()
#     # pool.join() 
            
#     # pool = Pool(n_jobs)
    
    
#     for index, row in subsub_df.iterrows():
#         if row.subject == 'P01':
#             applywarp = fsl.ApplyWarp()
#             applywarp.inputs.field_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-field.nii.gz'
#             applywarp.inputs.in_file = Template
#             applywarp.inputs.out_file = OutTemplatePath+'sub-' + row.subject + '_ses-' + row.session + '_TemplateRegistered.nii.gz'
#             applywarp.inputs.ref_file = row.path
#             applywarp.inputs.interp = 'nn'
#             applywarp.run()
#             #pool.apply_async(applywarp.run)
            
#     # pool.close()
#     # pool.join()
#     # pool = Pool(n_jobs)
    
    
#     for index, row in subsub_df.iterrows():
#         if row.subject == 'P01':
#             applywarp = fsl.ApplyWarp()
#             applywarp.inputs.field_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-field.nii.gz'
#             applywarp.inputs.in_file = Atlas
#             applywarp.inputs.out_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_AtlasRegistered.nii.gz'
#             applywarp.inputs.ref_file = row.path
#             applywarp.inputs.interp = 'nn'
#             #pool.apply_async(applywarp.run)
#             applywarp.run()
            
            
#     # pool.close()
#     # pool.join()
            
    
def Apply_FNIRT(DataPath, AtlasesPath, TemplatesPath, OutAtlasesPath, OutTemplatePath, n_jobs):

    layout = bids.BIDSLayout(DataPath, validate = False)
    layout_df = layout.to_df()
    sub_df = layout_df[(layout_df.extension == 'nii.gz')]
    subsub_df = sub_df[(sub_df.session == 'J0')]
    
    Template_layout = bids.BIDSLayout(TemplatesPath, validate = False)
    Template_layout_df = Template_layout.to_df()
    template_sub_df = Template_layout_df[(Template_layout_df.extension == 'nii.gz')]
    template_subsub_df = template_sub_df[(template_sub_df.session == 'J0')]
    
    Atlas_layout = bids.BIDSLayout(AtlasesPath, validate = False)
    Atlas_layout_df = Atlas_layout.to_df()
    atlas_sub_df = Atlas_layout_df[(Atlas_layout_df.extension == 'nii.gz')]
    atlas_subsub_df = atlas_sub_df[(atlas_sub_df.session == 'J0')]
        # pool = Pool(n_jobs)
    
    
    for index, row in subsub_df.iterrows():
        #if (row.subject in {'Patient01', 'Patient02'}) or (row.subject=='Patient05' and row.session == 'D03') or (row.subject=='Patient0' and row.session == 'D01'):
        if row.subject == 'P01':
            
            fnt = fsl.FNIRT()
            tmp = template_subsub_df[(template_subsub_df.subject==row.subject)]
            fnt.inputs.in_file = tmp.path[tmp.index[0]]
            fnt.inputs.ref_file = row.path
            fnt.inputs.field_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-field.nii.gz'
            fnt.inputs.warped_file = OutTemplatePath+'sub-' + row.subject + '_ses-' + row.session + '_MaskedAtlasRegistered.nii.gz'
            fnt.run()
            #pool.apply_async(fnt.run)
            
            
    # pool.close()
    # pool.join() 
            
    # pool = Pool(n_jobs)
    
    
    # for index, row in subsub_df.iterrows():
    #     if row.subject == 'P01':
    #         applywarp = fsl.ApplyWarp()
    #         tmp = template_subsub_df[(template_subsub_df.subject==row.subject)]
    #         applywarp.inputs.field_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-field.nii.gz'
    #         applywarp.inputs.in_file = tmp.path[tmp.index[0]]
    #         applywarp.inputs.out_file = OutTemplatePath+'sub-' + row.subject + '_ses-' + row.session + '_TemplateRegistered.nii.gz'
    #         applywarp.inputs.ref_file = row.path
    #         applywarp.inputs.interp = 'nn'
    #         applywarp.run()
    #         #pool.apply_async(applywarp.run)
            
    # pool.close()
    # pool.join()
    # pool = Pool(n_jobs)
    
    
    for index, row in subsub_df.iterrows():
        if row.subject == 'P01':
            applywarp = fsl.ApplyWarp()
            atl = atlas_subsub_df[(atlas_subsub_df.subject==row.subject)]
            applywarp.inputs.field_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_transform-field.nii.gz'
            applywarp.inputs.in_file = atl.path[atl.index[0]]
            applywarp.inputs.out_file = OutAtlasesPath+'sub-' + row.subject + '_ses-' + row.session + '_AtlasRegistered.nii.gz'
            applywarp.inputs.ref_file = row.path
            applywarp.inputs.interp = 'nn'
            #pool.apply_async(applywarp.run)
            applywarp.run()
            
            
    # pool.close()
    # pool.join() 
            