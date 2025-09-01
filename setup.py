__author__ = 'mahajrod'

import os
import sys
from setuptools import setup, find_packages

dependencies = ['scipy', 'numpy', 'pandas', 'matplotlib', 'venn',
                'biopython', 'xmltodict', 'statsmodels', "ete3", "bcbio-gff"]
if sys.version_info[1] < 7:  # check minor version of python3
    dependencies += ["importlib_resources"]
"""    
#  scipy numpy pandas matplotlib  biopython xmltodict bcbio-gff statsmodels pyparsing ete3 venn

if sys.version_info[0] == 3: # major version
    dependencies += ["ete3", 'venn', ]
    if sys.version_info[1] < 7: # check minor version of python3
        dependencies += ["importlib_resources"]

elif sys.version_info[0] == 2:
    dependencies += ["ete2"]
    dependencies += ["importlib_resources"]
else:
    raise ValueError("ERROR!!! Unsupported python version: %s" % str(sys.version_info[0]))
"""
setup(name='routoolpa',
      version='0.119',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=dependencies,
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),)
