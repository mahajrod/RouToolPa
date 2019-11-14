__author__ = 'mahajrod'

import os
import sys
from setuptools import setup, find_packages

dependencies = ['scipy', 'numpy', 'pandas', 'matplotlib', 'matplotlib-venn',
                'biopython', 'xmltodict', 'bcbio-gff', 'statsmodels']
if sys.version_info[0] == 3:
    dependencies += ["ete3"]
elif sys.version_info[0] == 2:
    dependencies += ["ete2"]
else:
    raise ValueError("ERROR!!! Unsupported python version: %s" % str(sys.version_info[0]))

setup(name='RouToolPa',
      version='0.83',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=dependencies,
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),)
