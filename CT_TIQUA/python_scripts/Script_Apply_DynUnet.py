#!/usr/bin/env python
# coding: utf-8


from monai.inferers import sliding_window_inference
from monai.transforms import (
    Compose,
    CropForegroundd,
    CopyItemsd,
    LoadImaged,
    Orientationd,
    ScaleIntensityRanged,
    Spacingd,
    AsDiscreted,
    Invertd,
    SaveImaged,
    EnsureChannelFirstd
)

from monai.networks.nets import DynUNet

from monai.data import (
    Dataset, DataLoader,
    decollate_batch,
)


import torch




def ApplyDynUnet(infile, model_path, outfolder, device):
    test_transforms = Compose(
        [
            LoadImaged(keys="image"),
            EnsureChannelFirstd(keys="image"),
            ScaleIntensityRanged(
                keys="image", 
                a_min=-15, 
                #a_min=-50,
                a_max=100,
                #a_max=150, 
                b_min=-1.0, 
                b_max=1.0, 
                clip=True),
            CropForegroundd(keys="image", source_key="image"),
            Orientationd(keys="image", axcodes="RAS"),
            Spacingd(
                keys="image",
                pixdim=(1, 1, 1),
                mode="bilinear",
            ),
        ]
    )

    post_transforms = Compose(
        [
            Invertd(
                keys="pred",
                transform=test_transforms,
                orig_keys="image",
                meta_keys="pred_meta_dict",
                orig_meta_keys="image_meta_dict",
                meta_key_postfix="meta_dict",
                nearest_interp=False,
                to_tensor=True,
            ),
            AsDiscreted(keys="pred", argmax=True) ,
            SaveImaged(keys="pred",
                       meta_keys="pred_meta_dict",
                       output_dir=outfolder,
                       output_postfix="seg",
                       separate_folder = False,
                       resample=False),
        ]
    )


    model =DynUNet(
            spatial_dims=3,
            in_channels=1,
            out_channels=8,
            kernel_size=[3, 3, 3, 3, 3, 3],
            strides=[1, 2, 2, 2, 2, [2, 2, 1]],
            upsample_kernel_size=[2, 2, 2, 2, [2, 2, 1]],
            norm_name="instance",
            deep_supervision=False,
            res_block=True,
        ).to(device)
    model.load_state_dict(torch.load(model_path))


    #test_dict_ds = [{"image":infile}]
    test_dict_ds = [{"image": image_name} for image_name in zip([infile])]
    test_ds = Dataset(data=test_dict_ds, transform=test_transforms)

    test_loader = DataLoader(test_ds, batch_size=1)# num_workers=4)


    #idx=1
    model.eval()
    with torch.no_grad():
        for test_data in test_loader:
            test_inputs = test_data["image"].to(device)
            roi_size = (96, 96, 96)
            sw_batch_size = 16
            test_data["pred"] = sliding_window_inference(test_inputs, roi_size, sw_batch_size, model, mode="gaussian", overlap=0.8)

            test_data = [post_transforms(i) for i in decollate_batch(test_data)]
            #test_output = from_engine(["pred"])(test_data)
            #print("idx number = ", idx, "/", len(test_loader))
            #idx+=1

#infile = '/Volumes/Mac_Data/SUMOONE_Miroir/2023_Ovalie/OVALIE_DATA/CT_nii_gz_Resampled/452751_452751-Crane_Sans_IV-Crane_Sans_IV_dcm---GI-02---2016-12-04000000.nii.gz'
#model_path = '/Volumes/Mac_Data/SUMOONE_Miroir/Scripts/CT-TIQUA/CT_TIQUA/data/Trained_DynUnet.pth'
#outfolder = '/Volumes/Mac_Data/SUMOONE_Miroir/'
#device = 'cpu'
#ApplyDynUnet(infile, model_path, outfolder, device)