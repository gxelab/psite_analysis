import sys
import re
import gzip
import click
from itertools import chain, islice

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO


def smart_open(path):
    """open plain text or gzipped file"""
    if path[-2:] == 'gz':
        return gzip.open(path, 'rt')
    else:
        return open(path, 'rt')


def strip_version(tx_name):
    """Strip transcript version"""
    return re.sub('\.\d+$', '', tx_name)


def read_fasta(path, ignore_version=False):
    """Construct a dict of input fasta"""
    fa = dict()
    with smart_open(path) as f:
        for record in SeqIO.parse(f, 'fasta'):
            if ignore_version:
                record.id = strip_version(record.id)
            fa[record.id] = str(record.seq)
    return fa


def read_txinfo(path, sep='auto'):
    """Read transcript infomation"""
    if sep == 'auto':
        return pd.read_csv(path, sep=None, engine='python')
    else:
        return pd.read_csv(path, sep=sep)


def chunk_iter(iterable, size=1024):
    """
    iterate by fixed block size
    
    Yield an iterator for each block. If the last chunk do not have enough
    items, all the remaining items is return as the last chunk.
    reference: reclosedev, 2012, https://stackoverflow.com/a/8998040/3926543
    """
    it = iter(iterable)
    while True:
        chunk = islice(it, size)
        try:
            fist_element = next(chunk)
        except StopIteration:
            return
        yield chain((fist_element,), chunk)


WC_PAIRING = str.maketrans('ACGTN', 'TGCAN')
def rev_comp(s):
    """reverse complementation DNA string"""
    return s[::-1].translate(WC_PAIRING)


@click.command(context_settings=dict(
    help_option_names=['-h', '--help'], show_default=True))
@click.argument('path_ref', type=click.STRING)
@click.argument('path_bam', type=click.STRING)
@click.option('-i', '--ignore_txversion', is_flag=True, default=False,
              help='ignore trasncript version in ".\d+" format')
@click.option('-n', '--nts', type=click.INT, default=1,
              help='fanking nucleotides to consider at each side')
@click.option('-o', '--output', type=click.File('wt'), default=sys.stdout)
def extract_features(path_ref, path_bam, output, ignore_txversion=False, nts=1):
    """
    Extract and save features for model training and testing

    \b
    path_ref: reference transcriptome (fasta) matching the bam
    path_bam: alignments of RPFs to reference transcriptome
    """
    ref = read_fasta(path_ref, ignore_version=ignore_txversion)
    with pysam.AlignmentFile(path_bam) as bam:
        for align in bam:
            if align.reference_name not in ref:
                continue
            # .reference_start: 0-based leftmost coordinate
            # .reference_end: reference_end points to one past the last aligned residue.
            # so that reference_end - reference_start = reference_length
            seq = ref[align.reference_name]
            if align.reference_start - nts >= 0:
                seq_left = seq[(align.reference_start - nts):(align.reference_start + nts)]
            else:
                seq_left = seq[:(align.reference_start + nts)].rjust(2*nts, '-')
            if align.reference_end + nts <= len(seq):
                seq_right = seq[(align.reference_end - nts):(align.reference_end + nts)]
            else:
                seq_right = seq[(align.reference_end - nts):].ljust(2*nts, '-')
            if align.is_reverse:
                flank = list(rev_comp(seq_left + seq_right))
            else:
                flank = list(seq_left + seq_right)
            out = [align.query_name, align.reference_name, str(align.reference_start + 1),
                   str(align.query_alignment_length)] + flank
            print('\t'.join(out), file=output)
    return


if __name__ == '__main__':
    extract_features()
