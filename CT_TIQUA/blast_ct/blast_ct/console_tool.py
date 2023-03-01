import sys
import json
import os
import argparse
import pandas as pd
import shutil
from .trainer.inference import ModelInference, ModelInferenceEnsemble
from .train import set_device
from .read_config import get_model, get_test_loader
from .nifti.savers import NiftiPatchSaver
import nibabel as nib

def path(string):
    if os.path.exists(string):
        return string
    else:
        sys.exit(f'File not found: {string}')


def console_tool():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', metavar='input', type=path, help='Path to input image.', required=True)
    parser.add_argument('--output', metavar='output', type=str, help='Path to output image.', required=True)
    parser.add_argument('--ensemble', help='Whether to use the ensemble (slower but more precise)', type=bool,
                        default=False)
    parser.add_argument('--device', help='GPU device index (int) or \'cpu\' (str)', default='cpu')

    parse_args, unknown = parser.parse_known_args()
    if not (parse_args.input[-7:] == '.nii.gz' or parse_args.input[-4:] == '.nii'):
        raise IOError('Input file must be of type .nii or .nii.gz')

    if not (parse_args.output[-7:] == '.nii.gz'):
        raise IOError('Output file must be of type .nii.gz')

    install_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(install_dir, 'data/config.json'), 'r') as f:
        config = json.load(f)

    device = set_device(parse_args.device)
    if device.type == 'cpu':
        config['test']['batch_size'] = 32
    job_dir = '/tmp/blast_ct'
    os.makedirs(job_dir, exist_ok=True)
    test_csv_path = os.path.join(job_dir, 'test.csv')
    pd.DataFrame(data=[['im_0', parse_args.input]], columns=['id', 'image']).to_csv(test_csv_path, index=False)

    model = get_model(config)
    test_loader = get_test_loader(config, model, test_csv_path, use_cuda=not device.type == 'cpu')
    saver = NiftiPatchSaver(job_dir, test_loader, write_prob_maps=False)

    if not parse_args.ensemble:
        model_path = os.path.join(install_dir, 'data/saved_models/model_1.pt')
        ModelInference(job_dir, device, model, saver, model_path, 'segmentation')(test_loader)
    else:
        model_paths = [os.path.join(install_dir, f'data/saved_models/model_{i:d}.pt') for i in range(1, 13)]
        ModelInferenceEnsemble(job_dir, device, model, saver, model_paths, task='segmentation')(test_loader)
    output_dataframe = pd.read_csv(os.path.join(job_dir, 'predictions/prediction.csv'))

    shutil.copyfile(output_dataframe.loc[0, 'prediction'], parse_args.output)
    shutil.rmtree(job_dir)

def console_tool_stand_alone(inp, out, device, prob_maps, ensemble, fold_tmp):
    #parser = argparse.ArgumentParser()
    #parser.add_argument('--input', metavar='input', type=path, help='Path to input image.', required=True)
    #parser.add_argument('--output', metavar='output', type=str, help='Path to output image.', required=True)
    #parser.add_argument('--ensemble', help='Whether to use the ensemble (slower but more precise)', type=bool,
    #                    default=False)
    #parser.add_argument('--device', help='GPU device index (int) or \'cpu\' (str)', default='cpu')

    #parse_args, unknown = parser.parse_known_args()
    #parse_args.input = inp
    #parse_args.output = out
    #parse_args.device = device
    #parse_args.ensemble = ensemble

    if not (inp[-7:] == '.nii.gz' or inp[-4:] == '.nii'):
        raise IOError('Input file must be of type .nii or .nii.gz')

    if not (out[-7:] == '.nii.gz'):
        raise IOError('Output file must be of type .nii.gz')

    install_dir = os.path.dirname(os.path.realpath(__file__))
    #with open(os.path.join(install_dir, 'data/config.json'), 'r') as f:
    with open(os.path.join(install_dir, 'data/config_7_classes.json'), 'r') as f:
        config = json.load(f)

    device = set_device(device)
    if device.type == 'cpu':
        config['test']['batch_size'] = 32
    #job_dir = '/tmp/blast_ct'
    job_dir = fold_tmp+'blast_ct'
    os.makedirs(job_dir, exist_ok=True)
    test_csv_path = os.path.join(job_dir, 'test.csv')
    pd.DataFrame(data=[['im_0', inp]], columns=['id', 'image']).to_csv(test_csv_path, index=False)

    model = get_model(config)
    test_loader = get_test_loader(config, model, test_csv_path, use_cuda=not device.type == 'cpu')
    saver = NiftiPatchSaver(job_dir, test_loader, write_prob_maps=True)

    if not ensemble:
        #model_path = os.path.join(install_dir, 'data/saved_models/model_1.pt')
        model_path = os.path.join(install_dir, 'data/7classes_models_pt/TL_REBLAST1_S5000.pt')
        ModelInference(job_dir, device, model, saver, model_path, 'segmentation')(test_loader)
    else:
        #model_paths = [os.path.join(install_dir, f'data/saved_models/model_{i:d}.pt') for i in range(1, 13)]
        model_paths = [os.path.join(install_dir, f'data/7classes_models_pt/TL_REBLAST{i:d}_S5000.pt') for i in range(1, 13)]
        ModelInferenceEnsemble(job_dir, device, model, saver, model_paths, task='segmentation')(test_loader)
    output_dataframe = pd.read_csv(os.path.join(job_dir, 'predictions/prediction.csv'))

    shutil.copyfile(output_dataframe.loc[0, 'prediction'], out) # déplacer la segmentation blast finale
    h = nib.load(output_dataframe.loc[0, 'prob_maps']) # le nifti contenant les 5 probability maps
    V = h.get_fdata()
    for ind in range(V.shape[4]):
        p_map = V[:,:,:,0,ind]
        if prob_maps[-4:] == '.nii':
            name = prob_maps[:-4]+str(ind)+prob_maps[-4:]
        elif prob_maps[-7:] == '.nii.gz':
            name = prob_maps[:-7]+str(ind)+prob_maps[-7:]
        else:
            raise IOError('ProbMap file must be of type .nii or .nii.gz')
        out_h = nib.Nifti1Image(p_map, h.affine)
        nib.save(out_h, name)
    #shutil.rmtree(job_dir)
    #shutil.rmtree(job_dir, ignore_errors=True)



