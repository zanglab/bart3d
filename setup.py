# Time-stamp: <2020-08-07>
'''
Copyright (c) 2020, 2021 Zhenjia Wang <zhenjia@virginia.edu>, Chongzhi Zang <zang@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: v1.0
@author: Zhenjia Wang, Chongzhi Zang
@contact: zhenjia@virginia.edu, zang@virginia.edu
'''

import sys
from setuptools import setup, find_packages
import bart3d

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")

# read requirements
with open('requirements.txt') as f:
    required = f.read().splitlines()

def main():
    if float(sys.version[:3])<3.0:
        sys.stderr.write("CRITICAL: Python version must be higher than or equal to 3.0!\n")
        sys.exit(1)
        
    setup(name="bart3d",
          version=bart3d.__version__,
          description="Inferring transcriptional regulators from differential Hi-C data",
          long_description=long_descr,
          author='Zhanjia Wang, Chongzhi Zang',
          author_email='zhenjia@virginia.edu, zang@virginia.edu',
          url='https://github.com/zanglab/bart3d',
          packages=find_packages(exclude=['testdata']),
          package_data={'':['bart3d.conf',
          				    'utility/bedGraphToBigWig',
          					'utility/*.sizes'],},
          include_package_data=True,
          scripts=['bin/bart3d',],
          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Science/Research',
              'License ::',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python',
              ],
          install_requires=required,
          )




if __name__=='__main__': 
    main()         
