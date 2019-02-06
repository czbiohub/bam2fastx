# -*- coding: utf-8 -*-

# Import modified 'os' module with LC_LANG set so click doesn't complain
from bam2fastx.os_utils import os  # noqa: F401
import sys

import click

from bam2fastx.fasta import fasta

settings = dict(help_option_names=['-h', '--help'])


@click.group(options_metavar='', subcommand_metavar='<command>',
             context_settings=settings)
def cli():
    """
    Aguamenti creates reflow batches for AWS jobs from experiment IDs
    """
    pass


cli.add_command(fasta, name='fasta')


if __name__ == "__main__":
    cli(sys.argv[1:])
