#!/usr/bin/env python
##################################################
# File Name: command.py
# Author: Rui
# mail: rdong1989@gmail.com
# Created Time: Tue 08 Nov 2022 09:35:37 AM EST
################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Usage: DISPbind (<command> | --help | --version)
'''

from docopt import docopt
import sys
from .version import __version__

__author__ = 'Rui Dong (rdong@mgh.harvard.edu)'

help_doc = '''

Usage: DISPbind <command> [options]

Command:
    align            Map fastq files bwa
    bam2bw           Generate bigwig file
    island           Characterize DisP islands from peak and bigwig files
'''


def main():
    # parse command
    command_log = 'DISPbind parameters: ' + ' '.join(sys.argv)
    if len(sys.argv) == 1:
        sys.exit(help_doc)
    elif sys.argv[1] == '--version' or sys.argv[1] == '-v':
        sys.exit(__version__)
    elif sys.argv[1] == 'align':
        from ./ import align
        align.align(docopt(align.__doc__, version=__version__))
    elif sys.argv[1] == 'bam2bw':
        from ./ import bam2bw
        bam2bw.bam2bw(docopt(bam2bw.__doc__, version=__version__))
    elif sys.argv[1] == 'island':
        from ./ import island
        island.island(docopt(island.__doc__, version=__version__))
    else:
        sys.exit(help_doc)


if __name__ == '__main__':
    main()
