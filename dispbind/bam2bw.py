'''
Usage: DISPbind.py bam2bw [options] -b bam -n NAME -o OUT

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -m MQ --mquality=MQUALITY      Mapping quality. [default: 10]
    -g GSIZE --gsize=GSIZE         Genome size file.
    -n NAME --name=NAME            Output file name. [default: bwa_out]
    -b BAM --bam=BAM               Input BAM file.
    -o OUT --output=OUT            Output directory. [default: alignment]
'''

import sys
import os
import os.path
import pysam
import pybedtools
from .helper import logger, which

__author__ = [
    'Rui Dong (rdong@mgh.harvard.edu)'
]

__all__ = ['align']

#@logger
def bam2bw(options):
    # check index files
    if not options['--bam']:
        sys.exit('Alignment requires Bam file!')
    if not options['--gsize']:
        sys.exit('Alignment requires genome size file!')

    # check output directory
    out_dir = check_outdir(options['--output'])
    generate_bw(out_dir, options['--name'], options['--mquality'], options['--bam'], options['--gsize'])

def check_outdir(out_dir):
    '''
    Create directory if not exists
    '''
    print('Check output directory...')
    # create output directory if not exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    dir_path = os.path.abspath(out_dir)
    return dir_path


def generate_bw(out_dir, name, mquality, bam, gsize):
    '''
    Create BigWig file
    '''
    # selece high quality reads
    filter_bam = 'samtools view -b -F2308 -q '
    filter_bam += ' %s %s > %s/%s ' % (mquality, bam, out_dir, name + '.filter.bam')
    return_code = os.system(filter_bam) >> 8
    if return_code:
        sys.exit('Error: cannot filter bam file!')
    # sort bam
    sort_bam = 'samtools sort -T '
    sort_bam += ' %s/%s -o %s/%s %s/%s ' % (out_dir, name + '.sorted', out_dir, name + '.sorted.bam', out_dir, name + '.filter.bam')
    return_code = os.system(sort_bam) >> 8
    if return_code:
        sys.exit('Error: cannot sort bam file!')
    # remove duplicate
    rm_dup = 'samtools rmdup -s '
    rm_dup += ' %s/%s %s/%s ' % (out_dir, name + '.sorted.bam', out_dir, name + '.sorted.deduped.bam')
    return_code = os.system(rm_dup) >> 8
    if return_code:
        sys.exit('Error: cannot remove dup!')
    # bam to bedpe and fill gaps
    bam2bedpe = 'samtools view -b -f2 '
    bam2bedpe += ' %s/%s | bedtools bamtobed -bedpe 2>/dev/null | cut -f 1,2,6 |sort -k1,1 -k2,2n > %s/%s' % (out_dir, name + '.sorted.deduped.bam', out_dir, name + '.bed')
    return_code = os.system(bam2bedpe) >> 8
    if return_code:
        sys.exit('Error: cannot remove dup!')

    # create Bigwig file
    if which('bedGraphToBigWig') is not None:
        print('Create BigWig file...')
        map_bam_fname = '%s/%s' % (out_dir, name + '.sorted.deduped.bam')
        # index bam if not exist
        if not os.path.isfile(map_bam_fname + '.bai'):
            pysam.index(map_bam_fname)
        map_bam = pysam.AlignmentFile(map_bam_fname, 'rb')
        # scale to HPB
        mapped_reads = map_bam.mapped
        s = 10000000.0 / mapped_reads
        map_bed_fname = '%s/%s' % (out_dir, name + '.bed')
        map_bed = pybedtools.BedTool(map_bed_fname)
        bedgraph_fname = '%s/%s' % (out_dir,name + '.bg')
        with open(bedgraph_fname, 'w') as bedgraph_f:
            for line in map_bed.genome_coverage(bg=True,
                                                g=gsize,
                                                scale=s, split=True):
                value = str(int(float(line[3]) + 0.5))
                bedgraph_f.write('\t'.join(line[:3]) + '\t%s\n' % value)
        # sort bedgraph
        sort_bg = 'LC_COLLATE=C sort -k1,1 -k2,2n %s/%s > %s/%s ' % (out_dir,name + '.bg', out_dir,name + '.sorted.bg')
        return_code = os.system(sort_bg) >> 8
        if return_code:
            sys.exit('Error: cannot sort bg!')
        # bg 2 bigwig
        bedgraph_sname = '%s/%s' % (out_dir,name + '.sorted.bg')
        bigwig_fname = '%s/%s' % (out_dir,name + '.bw')
        return_code = os.system('bedGraphToBigWig %s %s %s' %
                                (bedgraph_sname, gsize,
                                 bigwig_fname)) >> 8
        if return_code:
            sys.exit('Error: cannot convert bedGraph to BigWig!')
        file2rm = '%s/%s,%s/%s,%s/%s,%s/%s,%s/%s' % (out_dir, name + '.filter.bam', out_dir, name + '.sorted.bam', \
                                                     out_dir, name + '.bg', out_dir, name + '.sorted.bg', out_dir, name + '.bed')
        for f in file2rm.strip().split(","):
                os.remove(f)
    else:
        print('Could not find bedGraphToBigWig, so skip this step!')
