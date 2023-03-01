#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 15:52:21 2022

@author: clement
"""


import nibabel as nib
import nibabel.processing
import socket
import os
import shutil

from nipype.interfaces import fsl
import ants
import torch


# When executing from the commandline (install with pip)
from .blast_ct.blast_ct.console_tool import console_tool_stand_alone
from .python_scripts.Volume_estimation import Single_Volume_Inference

# When executing this script (from spyder for example)
#from blast_ct.blast_ct.console_tool import console_tool_stand_alone
#from python_scripts.Volume_estimation import Single_Volume_Inference


def inference(infile, outfolder, ensemble, keep_tmp_files):
   
    print('Start of the pipeline...')
    print('Summary:')
    print('infile='+infile)
    print('outfolder='+outfolder)
    print('ensemble='+str(ensemble))
    print('keep_tmp_files='+str(keep_tmp_files))
    sep = os.sep
    basename = os.path.basename(infile).split('.')[0]
    tmp_fold = outfolder+sep+'tmp'+sep
    os.makedirs(tmp_fold, exist_ok=True)
    
    fold = sep.join(os.path.realpath(__file__).split(sep)[:-1])
    

    matlab_App_path = sep+'compiled_matlab_scripts'+sep+'App'+sep+'application'+sep+'run_SkullStrip.sh'
    matlab_runtime_path = sep+'compiled_matlab_scripts'+sep+'RunTime'+sep+'v910'
    print('matlab_App_path='+matlab_App_path)
    print('matlab_runtime_path='+matlab_runtime_path)
    
    #CHECK THAT INPUT IMAGE HAS A QFORMCODE EQUAL TO 1
    print('Start of the quality control 1...')
    opt={"Sform_code":'aligned', "Qform_code":'unknown'}
    img_h = nib.load(infile)
    Vol = img_h.get_fdata()
    img_cleaned = nib.Nifti1Image(Vol, img_h.affine)
    sform_code = opt['Sform_code']
    qform_code = opt['Qform_code']
    #for i in img_h.header.items():
   # 	if i[0]=='sform_code':
   # 		initial_sform_code = i[1]
   # 	elif i[O]=='qform_code':
   # 		initial_qform_code = i[1]		
    img_cleaned.set_sform(img_cleaned.get_sform(), code=sform_code)
    img_cleaned.set_qform(img_cleaned.get_qform(), code=qform_code)
    nib.save(img_cleaned, tmp_fold+basename+'_clean.nii.gz')
    print('End of the quality control 1')
    
    #RESAMPLING
    print('Start of the resampling...')
    im_h = nib.load(tmp_fold+basename+'_clean.nii.gz')
    order = 0
    pixdim=[1,1,1]
    Im_resampled = nibabel.processing.resample_to_output(im_h, pixdim, order = order)
    resampled_file = tmp_fold+basename+'_Resampled.nii'
    nib.save(Im_resampled, resampled_file)
    print('End of the resampling')
    
    
    # #CHECK THAT RESAMPLED IMAGE HAS A QFORMCODE EQUAL TO 1
    
    
    
    #BRAIN EXTRACTION
    print('Start of the brain extraction...')
    #matlab_runtime_path = fold+sep+'matlab_scripts'+sep+'RunTime'+sep+'v910'
    #matlab_App_path = fold+sep+'matlab_scripts'+sep+'App'+sep+'application'+sep+'run_SkullStrip.sh'
    #print(matlab_runtime_path)
    outimage = tmp_fold+basename+'_SkullStripped.nii'
    outROI = tmp_fold+basename+'_ROI.nii'
    cmdline = matlab_App_path+' ' + matlab_runtime_path + ' ' + resampled_file + ' ' + outimage + ' ' + outROI
    #print(cmdline)
    os.system(cmdline)
    print('End of the brain extraction')
    
    
    #CHECK THAT SKULL STRIPPED AND ROI HAVE A QFORMCODE EQUAL TO 1
    
    print('Start of the quality control 2...')
    opt={"Sform_code":'scanner', "Qform_code":'scanner'}
    img_h = nib.load(tmp_fold+basename+'_SkullStripped.nii.gz')
    sform_code = opt['Sform_code']
    qform_code = opt['Qform_code']
    img_h.set_sform(img_h.get_sform(), code=sform_code)
    img_h.set_qform(img_h.get_qform(), code=qform_code)
    nib.save(img_h, tmp_fold+basename+'_SkullStripped_clean.nii.gz')
    print('End of the quality control 2')
    
    #SEGMENTATION BLAST
    print('Start of the segmentation...')
    segfile = outfolder+sep+basename+'_seg.nii.gz'
    probfile = tmp_fold+sep+basename+'_prob.nii.gz'
    if torch.cuda.is_available() and torch.cuda.device_count()>0:
        device = torch.cuda.current_device()
        print('Segmentation will run on GPU: ID='+str(device)+', NAME: '+torch.cuda.get_device_name(device))
    else:
        device = 'cpu'
        print('Segmentation will run on CPU')
    console_tool_stand_alone(resampled_file, segfile, device, probfile, ensemble, tmp_fold)
    print('End of the segmentation')
    
    #CHECK THAT SEGMENTATION HAS A QFORMCODE EQUAL TO 1
    
    # REGISTRATION
    print('Start of the linear registration...')
    Atlas = fold+sep+'data'+sep+'Resliced_Registered_Labels_mod.nii.gz'
    Template = fold+sep+'data'+sep+'TEMPLATE_miplab-ncct_sym_brain.nii.gz'
    flt = fsl.FLIRT()

    flt.inputs.in_file = Template
    flt.inputs.reference = tmp_fold+basename+'_SkullStripped_clean.nii.gz'
    flt.inputs.out_file = tmp_fold+basename+'_Template_FLIRTRegistered.nii'
    flt.inputs.out_matrix_file = tmp_fold+basename+ '_FLIRTRegisteredTemplate_transform-matrix.mat'
    flt.inputs.dof = 7
    flt.inputs.bins = 256
    flt.inputs.cost_func = 'normcorr'
    flt.inputs.interp = 'nearestneighbour'
    flt.inputs.searchr_x = [-180, 180]
    flt.inputs.searchr_y = [-180, 180]
    flt.inputs.searchr_z = [-180, 180]
    flt.run()
    


    applyxfm = fsl.ApplyXFM()
    applyxfm.inputs.in_matrix_file = tmp_fold+basename+ '_FLIRTRegisteredTemplate_transform-matrix.mat'
    applyxfm.inputs.in_file = Atlas
    applyxfm.inputs.out_file = tmp_fold+basename+'_Altas_FLIRTRegistered.nii'
    applyxfm.inputs.reference = tmp_fold+basename+'_SkullStripped_clean.nii.gz'
    applyxfm.inputs.apply_xfm = True
    applyxfm.inputs.out_matrix_file = tmp_fold+basename+ '_FLIRTRegisteredAtlas_transform-matrix.mat'
    applyxfm.inputs.interp = 'nearestneighbour'
    applyxfm.run()
            
    
    print('End of the linear registration')
    
    
    #CHECK THAT REGISTERED TEMPLATE AND ATLAS HAVE A QFORMCODE EQUAL TO 1
    
    
    
    print('Start of the elastic registration...')
    img_fixed = ants.image_read(tmp_fold+basename+'_SkullStripped_clean.nii.gz')
    img_moving = ants.image_read(tmp_fold+basename+'_Template_FLIRTRegistered.nii')
    outprefix=tmp_fold+basename
    reg = ants.registration(img_fixed, img_moving, outprefix=outprefix, random_seed=42)
    reg['warpedmovout'].to_file(tmp_fold+basename+'_Template_ANTSRegistered.nii.gz')
    
    mytx = reg['fwdtransforms']
    im_to_embarque = ants.image_read(tmp_fold+basename+'_Altas_FLIRTRegistered.nii')
    embarqued_im = ants.apply_transforms(img_fixed, im_to_embarque, transformlist=mytx, interpolator='nearestNeighbor')
    embarqued_im.to_file(outfolder+sep+basename+'_Altas_ANTSRegistered.nii.gz')
    print('End of the elastic registration')
    

    #CHECK THAT REGISTERED TEMPLATE AND ATLAS HAVE A QFORMCODE EQUAL TO 1

    print('Start of the volume computation...')
    seg = outfolder+sep+basename+'_seg.nii.gz'
    atlas = outfolder+sep+basename+'_Altas_ANTSRegistered.nii.gz'
    Labels = fold+sep+'data'+sep+'Labels_With_0.csv'
    outcsv = outfolder+sep+basename+'_Volumes.csv'
    Single_Volume_Inference(atlas, seg, Labels, outcsv)
    
    print('End of the volume computation')
    
    if not keep_tmp_files:
        print('Removing of the temporary files...')
        #shutil.rmtree(tmp_fold+'blast_ct'+sep)
        shutil.rmtree(tmp_fold, ignore_errors=True)
    
    print('End of the pipeline')






