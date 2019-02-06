# Import modified 'os' module with LC_LANG set so click doesn't complain
from .os_utils import os  # noqa: F401

from collections import defaultdict

import click


DELIMITER = "X"
FASTA_PREFIX = "aligned_sequences"
CELL_BARCODE = 'CB'
UMI = 'UB'

BAM_FILENAME = 'possorted_genome_bam.bam'
BARCODES_TSV = 'barcodes.tsv'


def read_single_column(filename):
    """Read single-column barcodes.tsv and genes.tsv files from 10x"""
    with open(filename) as f:
        lines = set(line.strip() for line in f)
    return lines


def read_10x_folder(folder):
    """Get QC-pass barcodes, genes, and bam file from a 10x folder

    Parameters
    ----------
    folder : str
        Name of a 10x cellranger output folder containing
        'possorted_genome_bam.bam' and 'barcodes.tsv' files

    Returns
    -------
    barcodes : list
        List of QC-passing barcodes from 'barcodes.tsv'
    bam_file : bamnostic.AlignmentFile
        Iterator over possorted_genome_bam.bam file
    """
    import bamnostic as bs

    barcodes = read_single_column(os.path.join(folder, BARCODES_TSV))

    bam_file = bs.AlignmentFile(os.path.join(folder, BAM_FILENAME), mode='rb')

    return barcodes, bam_file


def _pass_alignment_qc(alignment, barcodes):
    """Assert high quality mapping, QC-passing barcode and UMI of alignment"""
    high_quality_mapping = alignment.mapq == 255
    good_barcode = CELL_BARCODE in alignment.tags and \
        alignment.get_tag(CELL_BARCODE) in barcodes
    good_umi = UMI in alignment.tags

    pass_qc = high_quality_mapping and good_barcode and good_umi
    return pass_qc


def _parse_barcode_renamer(barcodes, barcode_renamer):
    """

    :param barcodes:
    :param barcode_renamer:
    :return:
    """
    if barcode_renamer is not None:
        renamer = {}

        with open(barcode_renamer) as f:
            for line in f.readlines():
                barcode, renamed = line.split()
                assert barcode in barcodes
                renamer[barcode] = renamed
    else:
        renamer = dict(zip(barcodes, barcodes))
    return renamer


def barcode_iterator(bam, barcodes, barcode_renamer):
    """Yield a (barcode, list of str) pair for each QC-pass barcode"""
    bam_filtered = (x for x in bam if _pass_alignment_qc(x, barcodes))

    renamer = _parse_barcode_renamer(barcodes, barcode_renamer)

    # alignments only have a CELL_BARCODE tag if they past QC
    bam_sort_by_barcode = sorted(bam_filtered,
                                 key=lambda x: x.get_tag(CELL_BARCODE))

    previous_barcode = None
    barcode_alignments = []
    for alignment in bam_sort_by_barcode:
        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        barcode = alignment.get_tag(CELL_BARCODE)

        # If this is a new non-null barcode, return all previous sequences
        if previous_barcode is not None and barcode != previous_barcode:
            yield renamer[previous_barcode], barcode_alignments

            # Reset the barcode alignments
            barcode_alignments = []

        # Add only the aligned sequence to this list of barcode alignments
        barcode_alignments.append(alignment.seq)

        # Set this current barcode as the previous one
        previous_barcode = barcode

    # Yield the final one
    yield renamer[previous_barcode], barcode_alignments


def _write_all_cells_in_one_file(cell_sequences, output_folder, fasta_prefix):
    filename = os.path.join(output_folder,
                            f"{fasta_prefix}.fasta")

    with open(filename, "w") as f:
        for cell, seq in cell_sequences.items():
            f.write(f">{cell}\n{seq}")
            # this "pass" makes PyCharm happy
            pass
    return filename


def _write_one_cell_per_file(cell_sequences, output_folder, fasta_prefix):
    os.makedirs(output_folder, exist_ok=True)

    filenames = []

    for cell, seq in cell_sequences.items():
        filename = os.path.join(output_folder, f"{fasta_prefix}_{cell}.fasta")
        with open(filename, "w") as f:
            f.write(f">{cell}\n{seq}")
        filenames.append(filename)
    return filenames


def write_cell_sequences(cell_sequences, output_folder,
                         one_cell_per_file=False, fasta_prefix=FASTA_PREFIX):
    if one_cell_per_file:
        filenames = _write_one_cell_per_file(cell_sequences, output_folder,
                                             fasta_prefix)
    else:
        filename = _write_all_cells_in_one_file(cell_sequences, output_folder,
                                                fasta_prefix)
        filenames = [filename]
    return filenames


def bam_to_fasta(bam, barcodes, barcode_renamer, output_folder, delimiter="X",
                 one_cell_per_file=False, fasta_prefix=FASTA_PREFIX):
    """Convert 10x bam to one-record-per-cell fasta

    Parameters
    ----------
    bam : bamnostic.AlignmentFile

    barcodes : list of str
        QC-passing barcodes
    barcode_renamer : str or None
        Tab-separated filename mapping a barcode to a new name, e.g.
        AAATGCCCAAACTGCT-1    lung_epithelial_cell|AAATGCCCAAACTGCT-1
    delimiter : str, default "X"
        Non-DNA or protein alphabet character to be ignored

    Returns
    -------
    filenames : list
        List of fasta filenames written
    """

    bam_filtered = (x for x in bam if _pass_alignment_qc(x, barcodes))

    renamer = _parse_barcode_renamer(barcodes, barcode_renamer)

    cell_sequences = defaultdict(str)

    for alignment in bam_filtered:

        # Get barcode of alignment, looks like "AAATGCCCAAACTGCT-1"
        barcode = alignment.get_tag(CELL_BARCODE)
        renamed = renamer[barcode]

        # Make a long string of all the cell sequences, separated
        # by a non-alphabet letter
        cell_sequences[renamed] += alignment.seq + delimiter + "\n"

    return write_cell_sequences(cell_sequences, output_folder,
                                one_cell_per_file, fasta_prefix)


@click.command()
@click.argument("tenx_folder")
@click.option('--all-cells-in-one-file/--one-cell-per-file', default=True,
              help="Create a single fasta, with each cell as a separate "
                   "record, whose sequences are separated by the delimiter "
                   f"'{DELIMITER}' (default), or create many fasta files, "
                   "one per cell")
@click.option('--barcode-renamer',
              help="Tab-separated file mapping barcodes (column 1) to renamed "
                   "ids (column 2)")
@click.option("--output-folder", help="Folder to output to. Default is "
                                      "current directory", default=".")
@click.option('--fasta-prefix', help="Filename prefix to use ",
              default=FASTA_PREFIX)
@click.option('--delimiter', default=DELIMITER)
def fasta(tenx_folder, all_cells_in_one_file, barcode_renamer=None,
          output_folder=".", fasta_prefix=FASTA_PREFIX, delimiter=DELIMITER):
    """Convert 10x bam to fasta of aligned sequences

    Parameters
    ----------
    tenx_folder : str
        Location of tenx folder containing possorted_genome_bam.bam and
        barcodes.tsv files
    """

    barcodes, bam = read_10x_folder(tenx_folder)
    one_cell_per_file = not all_cells_in_one_file

    filenames = bam_to_fasta(bam, barcodes, barcode_renamer=barcode_renamer,
                             output_folder=output_folder, delimiter=delimiter,
                             fasta_prefix=fasta_prefix,
                             one_cell_per_file=one_cell_per_file)
    if len(filenames) == 1:
        filename = filenames[0]
        click.echo(f"Wrote {filename}")
    else:
        n_files = len(filenames)
        click.echo(f"Wrote {n_files} fasta files in {output_folder}")
