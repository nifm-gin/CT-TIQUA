#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 17:00:19 2021

@author: cbrossard
"""


import bids
import nibabel as nib
import csv
import numpy as np
import nibabel.processing

def Computation_Volume_Seg_Atlas(Path_atlas, Path_Seg, label_file, Output_file_path, nb_classes_seg):
    
    layout_seg = bids.BIDSLayout(Path_Seg, validate = False)
    layout_seg_df = layout_seg.to_df()
    sub_seg_df = layout_seg_df[(layout_seg_df.extension == 'nii.gz')]
    
    layout_atl = bids.BIDSLayout(Path_atlas, validate = False)

    with open(label_file, newline='') as f:
        reader = csv.reader(f)
        Names_labels = list(reader)

    Lines = []
    
    for index, row in sub_seg_df.iterrows():
        if row.path == sub_seg_df.path[sub_seg_df.index[0]]: #If it is the first iteration
            Line_def = []
            Line_def.append("Subject")
            Line_def.append("Session")
            Line_def.append("Atlas")
            Line_def.append("Segmentation")
        seg_file = row.path
        atlas_file = layout_atl.get(subject = row.subject, session= 'J0', return_type='filename', extension='nii.gz')
        
        seg_h = nib.load(seg_file)
        atlas_h = nib.load(atlas_file[0])
        
        # The following if statement is abad patch to assure that the atlas and segmentation have the same shape.
        # The problem is in the resampling of the registered raw data at 1mm3 and then should be fixed on that function, not here ...
        if seg_h.shape != atlas_h.shape:
            print("Error")
            print(row)
            seg_h = nibabel.processing.conform(seg_h,atlas_h.shape, atlas_h.header.get("pixdim")[1:4], order = 0)
        else:
            print(row)
        
        seg = seg_h.get_fdata()
        atlas = atlas_h.get_fdata()
        
        
        Line = []
        Line.append(row.subject)
        Line.append(row.session)
        Line.append(Path_atlas)
        Line.append(Path_Seg)
        
        labels_values = np.unique(atlas)
        #seg_values = np.unique(seg)
        seg_values = range(nb_classes_seg)
        for ind, seg_val in enumerate(seg_values):
            mask_seg = seg==seg_val

            for i, label in enumerate(labels_values):
                if row.path == sub_seg_df.path[sub_seg_df.index[0]]: #If it is the first iteration
                    name = "Volume_seg"+str(seg_val)+"_"+ Names_labels[i][1]
                    Line_def.append(name)
                
                mask_atl = atlas==label

                intersect = mask_seg & mask_atl
                nb_vox = np.sum(intersect)
                Line.append(nb_vox)
                
        if row.path == sub_seg_df.path[sub_seg_df.index[0]]: #If it is the first iteration
            Lines.append(Line_def)
            
        Lines.append(Line)
        
    with open(Output_file_path, 'w', newline='') as csvfile:
        wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        for l in Lines:
            wr.writerow(l)
        
        
        
        
def compute_metrics(file_in, file_out):
    nb_regions_atlas = 19
    flag = 1
    with open(file_in, 'r', newline='') as infile:
        I = csv.reader(infile)
        with open(file_out, 'w', newline='') as outfile:
            O = csv.writer(outfile)
            for l, row in enumerate(I):
                if l==0:
                    Column_names = row
                    DataColumn_names = row[4:]
                    continue
                MetaData = row[:4]
                Data = row[4:]
                nb_seg = int(len(Data)/nb_regions_atlas)
                matrice = np.zeros((nb_seg, nb_regions_atlas))
                Names_types = []
                for s in range(nb_seg):
                    matrice[s,:] = Data[int(s*nb_regions_atlas):int((s+1)*nb_regions_atlas)]
                    Names_types.append('SUMMARY_seg'+str(s))
                Metrics_type = np.sum(matrice, axis=1)
                Metrics_Vol_atlas = np.sum(matrice, axis=0)
                Metrics_loc_vol = np.sum(matrice[1:,:], axis=0)
                Metrics_loc_prop = Metrics_loc_vol/Metrics_Vol_atlas
                matrice_proportion = np.zeros((nb_seg, nb_regions_atlas))
                for s in range(nb_seg):
                    matrice_proportion[s,:] = Data[int(s*nb_regions_atlas):int((s+1)*nb_regions_atlas)]
                    matrice_proportion[s,:] = matrice_proportion[s,:]/Metrics_Vol_atlas
                Names_full_proportion = [n+'_prop' for n in DataColumn_names]
                Full_data_prop = np.reshape(matrice_proportion, len(Names_full_proportion))
                Names_Vol_atlas = []
                Names_Prop_atlas = []
                for n in DataColumn_names[:nb_regions_atlas]:
                    splt = n.split('Volume_seg0_')[-1]
                    Names_Vol_atlas.append('SUMMARY_'+splt)
                    Names_Prop_atlas.append('Proportion_'+splt)
                
                
                OUTPUT_Column_names = Column_names+Names_types+Names_Vol_atlas+Names_Prop_atlas+Names_full_proportion
                OUTPUT_DATA = MetaData+Data+list(Metrics_type)+list(Metrics_loc_vol)+list(Metrics_loc_prop)+list(Full_data_prop)
                
                if flag:
                    O.writerow(OUTPUT_Column_names)
                    flag=0
                
                O.writerow(OUTPUT_DATA)
                
        #     Data.append(row)
        # wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        # for l in Lines:
        #     wr.writerow(l)
    
    
    
def Single_Volume_Inference(atlas, seg, Labels, outcsv):

    

    Names_lesions=['No_lesion','IPH', 'SDH', 'EDH', 'IVH', 'SAH', 'Petechiae', 'Edema']
    with open(Labels, newline='') as f:
        reader = csv.reader(f)
        Names_labels = list(reader)
        
    Zones_Numbers = list(range(len(Names_labels)))
    seg_h = nib.load(seg)
    atlas_h = nib.load(atlas)
        
    # The following if statement is abad patch to assure that the atlas and segmentation have the same shape.
    # The problem is in the resampling of the registered raw data at 1mm3 and then should be fixed on that function, not here ...
    if seg_h.shape != atlas_h.shape:
        print("Error")
        print("ATLAS SHAPE")
        print(atlas_h.shape)
        print("SEG SHAPE")
        print(seg_h.shape)
        print("PATH")
        print(atlas)
        seg_h = nibabel.processing.conform(from_img=seg_h,out_shape=atlas_h.shape, voxel_size=atlas_h.header.get("pixdim")[1:4], order = 0)
        
    seg = seg_h.get_fdata()
    atlas = atlas_h.get_fdata()
        

    Lines = []

    
    labels_values = np.unique(atlas)
    seg_values = range(8)
    for ind, seg_val in enumerate(seg_values):
        mask_seg = seg==seg_val
        name = "Volume_seg"+str(seg_val)
        Line = []
        Line.append(name)
        Line.append(Names_lesions[ind])
        for i, label in enumerate(labels_values):
            
            
            mask_atl = atlas==label

            intersect = mask_seg & mask_atl
            nb_vox = np.sum(intersect)
            Line.append(nb_vox)
        
        Lines.append(Line)


    with open(outcsv, 'w', newline='') as csvfile:
        wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        Just_label_names = [lab[1] for lab in Names_labels]
        Just_label_number = ['Volume_zone'+str(lab[0]) for lab in Names_labels]
        wr.writerow(['', '']+Just_label_number)
        wr.writerow(['', 'Name_Lesion']+Just_label_names)
        for l in Lines:
            wr.writerow(l)


   
def combined_regions_lesion_nifti(atlas, seg, atlas_labels_csv, out_nifti_path):

    seg_h = nib.load(seg)

    #print("segmentation nifti header: ")

    atlas_h = nib.load(atlas)
        
    # The following if statement is abad patch to assure that the atlas and segmentation have the same shape.
    # The problem is in the resampling of the registered raw data at 1mm3 and then should be fixed on that function, not here ...
    if seg_h.shape != atlas_h.shape:
        #print("Error between seg and atlas shapes")

        seg_h = nibabel.processing.conform(from_img=seg_h,out_shape=atlas_h.shape, voxel_size=atlas_h.header.get("pixdim")[1:4], order = 0)
        
    seg_data = seg_h.get_fdata()
    atlas_data = atlas_h.get_fdata()

    # Original atlas regions
    atlas_labels = []
    with open(atlas_labels_csv, newline='') as f:
        reader = csv.reader(f)
        atlas_labels = list(reader)
    #print(atlas_labels)
    lesion_names = ['No_lesion','IPH', 'SDH', 'EDH', 'IVH', 'SAH', 'Petechiae', 'Edema']

    # Combined atlas regions (supraregions)
    supratentorial = ["Right_frontal_cortex","Left_frontal_cortex","Right_temporal_cortex","Left_temporal_cortex","Right_insular_cortex","Left_insular_cortex",
                  "Right_parietal_cortex","Left_parietal_cortex","Right_occipital_cortex","Left_occipital_cortex","Right_basal_ganglia","Left_central_ganglia",
                  "Right_hypothalamus","Left_hypothalamus","Right_frontal_white_matter","Left_frontal_white_matter","Right_temporal_white_matter",
                  "Left_temporal_white_matter","Right_occipital_white_matter","Left_occipital_white_matter","Right_parietal_white_matter","Left_parietal_white_matter",
                  "Right_insular_white_matter","Left_insular_white_matter"]

    brainstem = ["Brainstem"]

    infratentorial = ["Right_cerebellum_lobe", "Left_cerebellum_lobe","Cerebellar_vermis"]

    extracerebral = ["Right_posterior_inferior_extra_cerebral_space","Left_posterior_inferior_extra_cerebral_space","Right_anterior_inferior_extra_cerebral_space",
                     "Left_anterior_inferior_extra_cerebral_space","Right_posterior_superior_extra_cerebral_space","Left_posterior_superior_extra_cerebral_space",
                     "Right_Superior_Anterior_Extra_Brain_Space","Left_Superior_Anterior_Extra_Brain_Space", "Rest"]


    final_nifti_labels = ['supratentorial_IPH', 'supratentorial_SAH', 'supratentorial_Petechiae', 'supratentorial_Edema',
                           'infratentorial_IPH', 'infratentorial_SAH', 'infratentorial_Petechiae', 'infratentorial_Edema',
                           'brainstem_IPH', 'brainstem_SAH', 'brainstem_Petechiae', 'brainstem_Edema', 
                           'SDH', 'EDH']
    

    new_combined_seg_nifti_data = np.zeros_like(seg_data)

    #print(np.unique(seg_data))

    for i in range(seg_data.shape[0]):
        for j in range(seg_data.shape[1]):
            for k in range(seg_data.shape[2]):

                if int(seg_data[i,j,k]) != 0: #0 is when there is no lesion
                    atlas_region_name = atlas_labels[int(atlas_data[i,j,k])][1]
                    lesion_name = lesion_names[int(seg_data[i,j,k])]
                    #print(f"{atlas_region_name}: {lesion_names[int(seg_data[i,j,k])]}")

                    combined_region_name = "-"
                    if atlas_region_name in supratentorial:
                        combined_region_name = "supratentorial"
                    elif atlas_region_name in brainstem:
                        combined_region_name = "brainstem"
                    elif atlas_region_name in infratentorial:
                        combined_region_name = "infratentorial"
                    elif atlas_region_name in extracerebral:
                        combined_region_name = ""
                    
                    if lesion_name  in ['SDH', 'EDH']: # cas a part pour ces lesions car elles sont pas rattache a des regions
                            final_nifti_label = lesion_name
                            final_nifti_index = final_nifti_labels.index(final_nifti_label)
                            new_combined_seg_nifti_data[i,j,k] = final_nifti_index
                    elif combined_region_name != "-":
                        

                        if combined_region_name != "":
                            final_nifti_label = combined_region_name+"_"+lesion_name
                        else:
                            final_nifti_label = lesion_name

                        #print(final_nifti_label)
                        if final_nifti_label in final_nifti_labels:
                            
                            final_nifti_index = final_nifti_labels.index(final_nifti_label)
                        
                            new_combined_seg_nifti_data[i,j,k] = final_nifti_index+1 #+1 because 0 corresponds to no label

    new_combined_seg_nifti = nib.Nifti1Image(new_combined_seg_nifti_data, seg_h.affine, seg_h.header)
    #print(np.unique(new_combined_seg_nifti_data))
    nib.save(new_combined_seg_nifti, out_nifti_path)
        


def ventricle_volume_computation(seg_path, out_csv_path):

    ventricle_seg = nib.load(seg_path)
    ventricle_seg_data = ventricle_seg.get_fdata()

    ventricle_volume = 0

    for i in range(ventricle_seg_data.shape[0]):
        for j in range(ventricle_seg_data.shape[1]):
            for k in range(ventricle_seg_data.shape[2]):
                    ventricle_volume += ventricle_seg_data[i,j,k]

    with open(out_csv_path, 'w', newline='') as csvfile:
        wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        wr.writerow(["ventricles_volume"])
        wr.writerow([ventricle_volume])
