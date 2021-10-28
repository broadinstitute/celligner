from setuptools import setup, find_packages
import sys
import os
import io
import subprocess

if sys.version_info.major < 3 or sys.version_info.minor < 2:
  raise ValueError("celligner is only compatible with Python 3.3 and above")
if sys.version_info.minor < 5:
  import warnings
  warnings.warn("celligner may not function properly on Python < 3.5")

#os.system('git submodule init && git submodule sync')

print("trying to install the required limma R package")
try:
  subprocess.run(
    'R -e \'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos="http://cran.us.r-project.org")};BiocManager::install("limma");\'', shell=True, check=True, 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except:
  print('failed to install limma. \
    please install R or check your R installation and then install limma with:\
    R -e \"if(!requireNamespace(\"BiocManager\", quietly = TRUE)){\
        install.packages(\"BiocManager\", repos=\"http://cran.us.r-project.org\")};\
      BiocManager::install(c(\"limma\"));\"')

print("Finished!")
def read(*paths, **kwargs):
  """Read the contents of a text file safely.
  >>> read("celligner", "VERSION")
  '0.1.0'
  >>> read("README.md")
  ...
  """

  content = ""
  with io.open(
    os.path.join(os.path.dirname(__file__), *paths),
    encoding=kwargs.get("encoding", "utf8"),
  ) as open_file:
    content = open_file.read().strip()
  return content


def read_requirements(path):
  return [
    line.strip()
    for line in read(path).split("\n")
    if not line.startswith(('"', "#", "-", "git+"))
  ]


setup(
  name='celligner',
  version=read("celligner", "VERSION"),
  description='A useful module for alligning cell lines to tumors',
  long_description=read("README.md"),
  long_description_content_type="text/markdown",
  author='Jeremie Kalfon',
  author_email='jkobject@gmail.com',
  url="https://github.com/BroadInstitute/celligner",
  packages=find_packages(exclude=["tests", ".github"]),
  package_data={'celligner': ['data/*']},
  python_requires='>=3.5',
  install_requires=read_requirements("requirements.txt"),
  entry_points={
    "console_scripts": ["celligner = celligner.__main__:main"]
  },
  extras_require={"test": read_requirements("requirements-test.txt")},
  classifiers=[
    "Programming Language :: Python :: 3",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
  ],
)

try: 
  subprocess.run(
    "pip install git+https://github.com/jkobject/mnnpy", shell=True, check=True, 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except:
  print('failed to install mnnpy. \
    please install Python or check your Python installation and then install mnnpy with:\
    pip install git+https://github.com/jkobject/mnnpy')