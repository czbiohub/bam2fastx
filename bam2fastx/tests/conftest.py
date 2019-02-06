import os

import pytest


@pytest.fixture
def data_folder():
    """Absolute path to where test data is stored"""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        './data')


@pytest.fixture
def tenx_folder(data_folder):
    """Absolute path to where test 10x bam file and barcodes file are stored"""
    return os.path.join(data_folder, "10x-example")
