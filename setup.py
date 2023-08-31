from setuptools import setup, find_packages
from codecs import open
from os import path

this_dir = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(this_dir, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


LONG_DESCRIPTION = """
**autoGDC** is a Python package that provides a simple way to interface with
the wonderful Genomic Data Commons repository.
"""

setup(
    name='autoGDC',
    version='0.0.1',
    description='Automatic Genomic Data Commons processing for bioinformaticians',
    long_description=LONG_DESCRIPTION,
    url='https://gitlab.com/chase_brown/autogdc',
    author='Chase Brown',
    author_email='chase.brown.2016@gmail.com',

    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 1 - Planning',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # License employed
        'License :: OSI Approved :: MIT License',

        # Python versions supported
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics meta-analysis GDC Genomic-Data-Commons data-repository',

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires=[
                      "pandas>=1.2", #df.to_parquet() needs to work w/o path
                      "diskcache",
                      "prefixed",
                      "tables",
                      "colorama",
                      "termcolor", # for gdc-client
                      "intervaltree", # for gdc-client
                      "rpy2",
                      "mygene",
                      "requests_cache",
                      "geode @ git+git://github.com/Maayanlab/geode.git"
                      ]
)
