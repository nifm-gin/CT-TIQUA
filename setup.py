import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CT_TIQUA",
    version="2.0",
    author="Brossard Clement, Lemasson Benjamin",
    author_email="benjamin.lemasson@univ-grenoble-alpes.fr",
    description="computed Tomography traumatic brain Injury QUAntification",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/brosscle/CT-TIQUA",
    packages=['CT_TIQUA', 'CT_TIQUA.blast_ct', 'CT_TIQUA.blast_ct.blast_ct', 'CT_TIQUA.blast_ct.blast_ct.trainer', 'CT_TIQUA.blast_ct.blast_ct.nifti', 'CT_TIQUA.blast_ct.blast_ct.models', 'CT_TIQUA.python_scripts'],
    package_data={'': ['data/*', 'README.md', 'data/Labels_With_0.csv', 'data/Resliced_Registered_Labels_mod.nii.gz', 'data/TEMPLATE_miplab-ncct_sym_brain.nii.gz', 'data/7classes_models_pt/TL_REBLAST*_S5000.pt']},
    scripts=['CT_TIQUA/python_scripts/Volume_estimation.py'],
    entry_points={
        'console_scripts': [
            'ct-tiqua = CT_TIQUA.console_tool:console_tool',
        ]
    },
    install_requires=[
        'scipy>=1.4.0',
        'numpy>=1.21',
        'pandas',
        'nibabel',
        'torch',
        'SimpleITK==1.2.4',
        'tensorboard',
        'pybids',
        'antspyx',
        'nipype',
        'boutiques',
    ],
    python_requires='>=3.7',
)
