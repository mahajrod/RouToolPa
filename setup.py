__author__ = 'mahajrod'
import os
from setuptools import setup, find_packages

dependencies = ['scipy', 'numpy', 'pandas', 'matplotlib', 'matplotlib-venn',
                'biopython', 'xmltodict', 'bcbio-gff', 'statsmodels']

setup(name='MAVR',
      version='0.80',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=dependencies,
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),)
