'''
Usage: DISPbind align [options] -i INDEX -a FQ1 -b FQ2 -o OUT -n NAME

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -i INDEX --index=INDEX         Index files for BWA
    -p THREAD --thread=THREAD      Running threads. [default: 10]
    -m MQ --mquality=MQUALITY      Mapping quality. [default: 10]
    -g GSIZE --gsize=GSIZE         Genome size file.
    -n NAME --name=NAME            Output file name. [default: bwa_out]
    -a FQ1 --fastq1=FQ1            Input R1 file.
    -b FQ2 --fastq2=FQ2            Input R2 file.
    -o OUT --output=OUT            Output directory. [default: alignment]
'''

import sys
import os
import os.path
import numpy.core.multiarray
import pysam
import pybedtools
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from helper import logger, which

__author__ = [
    'Rui Dong (rdong@mgh.harvard.edu)'
]

__all__ = ['align']

#@logger
def align(options):
    # check index files
    if not options['--index']:
        sys.exit('Alignment requires BWA index files!')
    if not options['--fastq1'] or not options['--fastq2']:
        sys.exit('Alignment requires paired-end fastq files!')
    if options['--fastq1'] == options['--fastq2']:
        sys.exit('Check the -a and -b input, they should not be the same file!')

    # check output directory
    out_dir = check_outdir(options['--output'])
    bwa_map(out_dir,options['--index'], options['--name'], options['--mquality'], options['--fastq1'], options['--fastq2'], options['--thread'], options['--gsize'])

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


def bwa_map(out_dir, index, name, mquality, fastq1, fastq2, thread, gsize):
    '''
    1. Map reads with BWA
    2. Create BigWig file
    '''
    print('Map reads with BWA...')
    bwa_cmd1 = 'bwa mem -t '
    bwa_cmd1 += ' %s %s %s %s' % (thread, index, fastq1, fastq2)
    bwa_cmd1 += '| samtools sort -T %s/%s -o %s/%s 2>/dev/null' % (out_dir, name + '.sort.tmp', out_dir, name + '.raw.bam')
    return_code = os.system(bwa_cmd1) >> 8
    if return_code:
        sys.exit('Error: cannot map reads with BWA!')

    filter_bam = 'samtools view -b -F2308 -q '
    filter_bam += ' %s %s/%s > %s/%s ' % (mquality, out_dir, name + '.raw.bam', out_dir, name + '.bam')
    return_code = os.system(filter_bam) >> 8
    if return_code:
        sys.exit('Error: cannot filter bam file!')
    # sort bam
    sort_bam = 'samtools sort -T '
    sort_bam += ' %s/%s -o %s/%s %s/%s ' % (out_dir, name + '.sorted.temp', out_dir, name + '.sorted.bam', out_dir, name + '.bam')
    return_code = os.system(sort_bam) >> 8
    if return_code:
        sys.exit('Error: cannot sort bam file!')
   # remove duplicate
    rm_dup = 'samtools rmdup -s '
    rm_dup += ' %s/%s %s/%s ' % (out_dir, name + '.sorted.bam', out_dir, name + '.deduped.bam')
    return_code = os.system(rm_dup) >> 8
    if return_code:
        sys.exit('Error: cannot remove dup!')
    # bam to sortn
    bamsortn = 'samtools view -f2 -b'
    bamsortn += ' %s/%s |samtools sort -T %s/%s -n > %s/%s' % (out_dir, name + '.deduped.bam',out_dir, name + '.sort.tmp', out_dir, name + '.sorted_n.bam')
    return_code = os.system(bamsortn) >> 8
    if return_code:
        sys.exit('Error: cannot sorted by name!')

    bam2bedpe = 'bedtools bamtobed -bedpe -i '
    bam2bedpe += ' %s/%s 2>/dev/null | perl -alne \'print \"$F[0]\\t$F[1]\\t$F[5]\" if $F[5]>$F[1] and $F[0]=~/chr/ and $F[1]>0 and $F[5]>0\' |sort -k1,1 -k2,2n > %s/%s' % (out_dir, name + '.sorted_n.bam', out_dir, name + '.bed')
    return_code = os.system(bam2bedpe) >> 8
    if return_code:
        sys.exit('Error: cannot remove dup!')


    # create Bigwig file
    if which('bedGraphToBigWig') is not None:
        print('Create BigWig file...')
        map_bam_fname = '%s/%s' % (out_dir, name + '.deduped.bam')
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
                value = str(float(line[3]) + 0.5)
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

        file2rm = '%s/%s,%s/%s,%s/%s,%s/%s,%s/%s' % (out_dir, name + '.bam',out_dir, name + '.sorted_n.bam', out_dir, name + '.bg', out_dir, name + '.sorted.bg', out_dir, name + '.bed')
        for f in file2rm.strip().split(","):
                os.remove(f)
    else:
        print('Could not find bedGraphToBigWig, so skip this step!')
