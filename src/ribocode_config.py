import sys
import re

if __name__:
    header = [
        '# SampleName', 'AlignmentFile', 'Stranded(yes/reverse)',
        'P-siteReadLength', 'P-siteLocations']
    print('\t'.join(header))

    alignment_file, config_file = sys.argv[1:]
    with open(config_file, 'rt') as fh:
        for line in fh:
            if line.startswith('# Read lengths used:'):
                qwidth_range = re.sub(r'.*?\[(.*?)\].*', r'\1', line.rstrip())
                break
    qwidth_range = re.sub(' ', '', qwidth_range)
    out = [
        alignment_file.removesuffix('.bam'), alignment_file,
        'yes', qwidth_range, re.sub(r'\d+', '12', qwidth_range)
    ]
    print('\t'.join(out))
