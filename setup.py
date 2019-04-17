__author__ = 'mahajrod'
import os
from setuptools import setup, find_packages

setup(name='MAVR',
      version='0.58',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=['scipy', 'numpy', 'pandas', 'matplotlib', 'matplotlib-venn', 'ete2', 'biopython', 'xmltodict', 'bcbio-gff', 'statsmodels'],
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),)
