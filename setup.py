from setuptools import setup
import sys
import os
if sys.version_info.major < 3 or sys.version_info.minor < 2:
    raise ValueError("Celligner is only compatible with Python 3.3 and above")
if sys.version_info.minor < 5:
    import warnings
    warnings.warn("Celligner may not function properly on Python < 3.5")

os.system('git submodule init && git submodule sync')

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='Celligner',
    version='1.0',
    description='A useful module for any CompBio',
    long_description=long_description,
    author='Jeremie Kalfon',
    author_email='jkobject@gmail.com',
    url="https://github.com/BroadInstitute/Celligner",
    packages=['celligner'],
    package_data={'celligner': ['data/*']},
    python_requires='>=3.5',
    install_requires=[
        #'genepy',
        ## from requirements.txt
        "numpy",
        "pandas",
        "scikit_learn",
        "umap_learn",
        "contrastive",
        "mnnpy",
        "snn",
        "umap",
    ],  # external packages as dependencies
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
