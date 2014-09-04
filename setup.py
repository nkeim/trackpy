#!/usr/bin/env python

from distutils.core import setup

setup(
    name='tasker',
    version='0.2',
    author='Nathan C. Keim',
    author_email='nkeim@seas.upenn.edu',
    url='https://github.com/nkeim/tasker',
    packages=['tasker'],
    install_requires=['path.py',],
    test_suite = 'nose.collector'
    )
