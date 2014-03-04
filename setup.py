from setuptools import setup, Extension, find_packages
import sys

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    raise SystemExit("PhasingTools requires Python 2.7")

globals = {}
execfile("src/pbhla/__init__.py", globals)
__VERSION__ = globals["__VERSION__"]

DESC = 'Tools for analyzing HLA data from SMRT sequencing'

required = [
    "numpy >= 1.7.0"
    "pbcore >= 0.8.0",
    "pypeflow >= 0.1.1",
    "pbtools.pbdagcon >= 0.2.3",
    "PhasingTools >= 0.1.0"
]

extras = {
    "ClassII": [
        "HBAR-DTK >= 0.1.5"
    ]
}

setup(
    name = 'HlaTools',
    version=__VERSION__,
    author='Brett Bowman',
    author_email='bbowman@pacificbiosciences.com',
    url='https://github.com/bnbowman/HlaTools',
    description=DESC,
    license=open('LICENSES.txt').read(),
    packages = find_packages('src'),
    package_dir = {'':'src'},
    zip_safe = False,
    install_requires = required,
    extras_require = extras
)
