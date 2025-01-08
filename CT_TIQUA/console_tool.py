#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 15:52:21 2022

@author: clement
"""

import sys
import os
import argparse
from .main import inference

def path(string):
    if os.path.exists(string):
        return string
    else:
        sys.exit(f'File not found: {string}')


def console_tool():
    #parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', metavar='input', type=path, help='Path to input image, as a nifti file, compressed of not (str).', required=True)
    parser.add_argument('--output', metavar='output', type=str, help='Path to the output folder where output data will be stored (str).', required=True)

    parser.add_argument('--keep_tmp_files', help='Keep temporary files at the end of the pipeline.', action="store_true", default=False)
    parser.add_argument('--ventricles', help='Enable ventricles volume comutation (TotalSegmentator).', action="store_true", default=False)

    
    parse_args, unknown = parser.parse_known_args()
    if not (parse_args.input[-7:] == '.nii.gz' or parse_args.input[-4:] == '.nii'):
        raise IOError('Input file must be of type .nii or .nii.gz')

    if (parse_args.output[-7:] == '.nii.gz' or parse_args.output[-4:] == '.nii'):
        raise IOError('Output must be a folder, not an image.')

    os.makedirs(parse_args.output+'/tmp/', exist_ok=True)
    inference(parse_args.input, parse_args.output, parse_args.keep_tmp_files, ventricles_seg=parse_args.ventricles)
