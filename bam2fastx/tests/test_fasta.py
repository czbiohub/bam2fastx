#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from click.testing import CliRunner
import pytest


@pytest.fixture
def fasta_folder(data_folder):
    return os.path.join(data_folder, "fasta")

@pytest.fixture
def all_cells_one_file(fasta_folder):
    from bam2fastx.fasta import FASTA_PREFIX
    return os.path.join(fasta_folder, "all-cells-in-one-file",
                        f"{FASTA_PREFIX}.fasta")

@pytest.fixture
def one_cell_per_file_folder(fasta_folder):
    return os.path.join(fasta_folder, "one-cell-per-file")

@pytest.fixture
def one_cell_per_file_listing(one_cell_per_file_folder):
    return sorted(os.listdir(one_cell_per_file_folder))


def test_fasta_all_cells_one_file(tenx_folder, all_cells_one_file, tmpdir):
    from bam2fastx.fasta import fasta, FASTA_PREFIX

    runner = CliRunner()
    result = runner.invoke(fasta, [tenx_folder, "--output-folder", tmpdir])

    # exit code of '0' means success!
    assert result.exit_code == 0
    assert 'Wrote' in result.output
    assert 'aligned_sequences.fasta' in result.output

    # Make sure the files are there
    fasta = os.path.join(tmpdir, FASTA_PREFIX + ".fasta")
    assert os.path.exists(fasta)

    with open(all_cells_one_file) as f:
        true_fasta = f.readlines()
    with open(fasta) as f:
        test_fasta = f.readlines()

    assert test_fasta == true_fasta



def test_fasta_one_cell_per_file(tenx_folder, one_cell_per_file_folder,
                                 one_cell_per_file_listing, tmpdir):
    from bam2fastx.fasta import fasta, FASTA_PREFIX

    runner = CliRunner()
    result = runner.invoke(fasta, [tenx_folder, "--output-folder", tmpdir,
                                   "--one-cell-per-file"])

    n_files = len(one_cell_per_file_listing)

    # exit code of '0' means success!
    assert result.exit_code == 0
    assert f'Wrote {n_files} fasta files' in result.output
    assert len(os.listdir(tmpdir)) == n_files

    for filename in one_cell_per_file_listing:
        f1 = os.path.join(one_cell_per_file_folder, filename)
        f2 = os.path.join(tmpdir, filename)

        # Ensure file contents are the same
        with open(f1) as f:
            true_fasta = f.readlines()
        with open(f2) as f:
            test_fasta = f.readlines()

        assert test_fasta == true_fasta

