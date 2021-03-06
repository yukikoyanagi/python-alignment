#!/usr/bin/env python
# coding=utf-8

from distutils.core import setup

# noinspection PyArgumentList
setup(
    name='alignment',
    version='1.1.11',
    author='Eser Aygün',
    author_email='eser.aygun@gmail.com',
    packages=['alignment'],
    url='https://github.com/eseraygun/python-alignment',
    license='BSD 3-Clause License',
    description='Native Python library for generic sequence alignment.',
    long_description=open('README.rst').read(),
    requires=['numpy', 'six'],
    package_dir={'alignment': 'alignment'},
    package_data={'alignment': ['config/clustermode.pkl',
                                'config/blosum62.txt',
                                'config/blosum80.txt',
                                'config/blosum90.txt',
                                'config/bstrands.txt']},
    scripts=['scripts/getscores.py',
             'scripts/comparewithgdt.py',
             'scripts/alignlocpats.py',
             'scripts/computediff.py']
)
