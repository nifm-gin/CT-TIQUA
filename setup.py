import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CT_TIQUA",
    version="4.0",
    author="Brossard Clement, Fehr Delude Theotime, Lemasson Benjamin",
    author_email="benjamin.lemasson@univ-grenoble-alpes.fr",
    description="Computed Tomography Traumatic brain Injury QUAntification",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nifm-gin/CT-TIQUA",
    packages=['CT_TIQUA', 'CT_TIQUA.python_scripts'],
    package_data={'': ['data/*', 'README.md', 'data/Labels_With_0.csv', 'data/Resliced_Registered_Labels_mod.nii.gz', 'data/TEMPLATE_miplab-ncct_sym_brain.nii.gz', 'data/ArterialAtlas.nii.gz', 'data/Labels_With_0_vasc.csv', 'data/24-02-23-13h18m_best_model.pt']},
    scripts=['CT_TIQUA/python_scripts/Volume_estimation.py', 'CT_TIQUA/python_scripts/Script_Apply_DynUnet.py'],
    entry_points={
        'console_scripts': [
            'ct-tiqua = CT_TIQUA.console_tool:console_tool',
        ]
    },
    install_requires=[
        'acvl-utils==0.2.2',
        'torch>=2.0.0',
        'totalsegmentator',
        'scipy>=1.4.0',
        'monai',
        'numpy>=1.21',
        'pandas',
        'nibabel',
        'SimpleITK>=1.2.4',
        'tensorboard',
        'pybids',
        'antspyx',
        'nipype',
        'boutiques',
    ],
    python_requires='>=3.10',
)
