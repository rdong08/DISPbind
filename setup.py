#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from dispbind.version import __version__

setup(name='DISPbind',
      version=__version__,
      description='Disorder protein genomic binding analysis toolkit',
      author='Rui Dong',
      author_email='rdong@mgh.harvard.edu',
      url='https://github.com/rdong08/DISPbind',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      keywords='disorder protein',
      packages=find_packages(),
      install_requires=[
          'requests',
          'pysam>=0.8.4',
          'pybedtools>=0.7.5',
          'pandas>=1.0.5',
          'matplotlib>=3.3.0',
          'numpy>=1.10.0',
          'docopt>=0.6.0',
          'pyBigWig>=0.3.17',
          'seaborn>=0.10.0',
          'scipy',
      ],
      entry_points={
          'console_scripts': [
              'DISPbind=dispbind.DISPbind:main'
          ],
      },
      )
