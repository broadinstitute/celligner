from setuptools import setup, find_packages
import sys
import os
import io

if sys.version_info.major < 3 or sys.version_info.minor < 2:
    raise ValueError("Celligner is only compatible with Python 3.3 and above")
if sys.version_info.minor < 5:
    import warnings
    warnings.warn("Celligner may not function properly on Python < 3.5")

os.system('git submodule init && git submodule sync')

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
    name='Celligner',
    version=read("Celligner", "VERSION"),
    description='A useful module for alligning cell lines to tumors',
    long_description=read("README.md"),
    author='Jeremie Kalfon',
    author_email='jkobject@gmail.com',
    url="https://github.com/BroadInstitute/Celligner",
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
