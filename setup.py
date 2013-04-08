from setuptools import setup, Extension, find_packages
import sys

desc = 'Tools for analyzing and phasing SMRT sequencing data from ' + \
       'large genonic amplicons, primarily HLA and MHC genes'

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    raise SystemExit("HlaTools requires Python 2.7")

setup(
    name = 'HlaTools',
    version='0.1.0',
    author='Brett Bowman',
    author_email='bbowman@pacificbiosciences.com',
    url='https://github.com/bnbowman/HlaTools',
    description=desc,
    license=open('LICENSES.txt').read(),
    packages = find_packages('src'),
    package_dir = {'':'src'},
    zip_safe = False,
    install_requires=[]
    )
