#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read()

test_requirements = [
    'pytest', 'coverage'
]

setup(
    name='bam2fastx',
    version='0.1.0',
    description="Convert 10x bam file to individual FAST{A,Q} files of aligned reads per cell",
    long_description=readme + '\n\n' + history,
    author="Olga Botvinnik",
    author_email='olga.botvinnik@czbiohub.org',
    url='https://github.com/olgabot/bam2fastx',
    packages=[
        'bam2fastx',
    ],
    package_dir={'bam2fastx':
                 'bam2fastx'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='bam2fastx',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    entry_points={
        'console_scripts': [
            'bam2fastx = bam2fastx.commandline:cli'
        ]
    },
    test_suite='tests',
    tests_require=test_requirements
)
