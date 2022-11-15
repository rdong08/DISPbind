'''
Usage: DISPbind.py island [options] -p PEAK (-b BW | -l LIST) -o OUT

Options:
    -h --help                             Show help message.
    -v --version                          Show version.
    -p PEAK --peak=PEAK                   DisP peaks.
    -b BW --bigwig=BW                     Bigwig file for corresponding sample.
    -o OUTPREFIX --output=OUTPREFIX       Output island file. [default: DisP]
    -l LIST --list=LIST                   Bigwig replicate (will calculate average signals for replicates)
    --plot                                Hockey plot for signal vs rank.
'''

import sys
import os
import os.path
import pyBigWig
import numpy as np
import pandas as pd
import argparse
import re
import os
from .helper import logger, which
from .signal import fetch_genome_size, fetch_bed_loc, fetch_signal
import seaborn as sns
import matplotlib.pyplot as plt


__author__ = [
    'Rui Dong (rdong@mgh.harvard.edu)'
]

__all__ = ['island']

#@logger
def island(options):
    # check index files
    if not options['--peak']:
        sys.exit('Island identification requires peak files!')
    if not options['--bigwig'] and not options['--list']:
        sys.exit('Island identification requires peak files!')

    # bedtools merge loops
    print('Sort and merge peaks of distance shorter than 20k...')
    sort_peaks = 'bedtools sort -i %s > %s' % (options['--peak'], options['--output'] + '.sort.peak')
    return_code = os.system(sort_peaks) >> 8
    if return_code:
        sys.exit('Error: cannot sort peaks!')
    merge_near_peaks = 'bedtools merge -d 20000 -i %s -c 4 -o collapse > %s' % (options['--output'] + '.sort.peak', options['--output'] + '.merged.peak')
    return_code = os.system(merge_near_peaks) >> 8
    if return_code:
        sys.exit('Error: cannot merge peaks!')
    merged_bed_file = '%s' % (options['--output'] + '.merged.peak')
    bed_loc = fetch_bed_loc(merged_bed_file)
    print('Start calling signals...')
    if options['--bigwig']:
        ##for each bigwig file, calculate signals of single replicate
        call_signal = fetch_signal(options['--bigwig'], bed_loc)
        df_island = pd.DataFrame(bed_loc, columns = ['chr','start','end'])
        df_island['signal'] = call_signal
    elif options['--list']:
        ##for each bigwig file, calculate signals of multiple replicate files
        bw_signal=[]
        sample_name=[]
        for line in open(options['--list']):
            cols = line.strip().split()
            if (len(cols) == 1):
                file_name = os.path.basename(cols[0])
                bw_out_name = os.path.splitext(file_name)[0]
                bw_file = cols[0]
            else:
                bw_out_name = cols[0]
                bw_file = cols[1]
            sample_name.append(bw_out_name)
            bw_signal_line = fetch_signal(bw_file, bed_loc)
            bw_signal.append(bw_signal_line)
        df_island = pd.DataFrame(bed_loc, columns = ['chr','start','end'])
        df_island['signal'] = np.average(bw_signal, axis = 0)

    else:
        sys.exit('Error: could not find bigwig files!')

    print('Identify DisP island...')
    df_island['rank'] = df_island['signal'].rank(method='max')
    #df_island['scale_signal'] = (df_island['signal'] - np.min(df_island['signal'])) / (np.max(df_island['signal']) - np.min(df_island['signal']))
    df_island['scale_signal'] = (df_island['signal']) / (np.max(df_island['signal']) - np.min(df_island['signal']))
    #df_island['scale_rank'] = (df_island['rank'] - np.min(df_island['rank'])) / (np.max(df_island['rank']) - np.min(df_island['rank']))
    df_island['scale_rank'] = (df_island['rank']) / (np.max(df_island['rank']) - np.min(df_island['rank']))
    df_island['diff'] = df_island['scale_rank'] - df_island['scale_signal']
    tangent = df_island.loc[df_island['diff'] == max(df_island['diff']),"signal"]
    df_island['DisP_island'] = 'non_island'
    df_island.loc[df_island['signal'] >= tangent.values[0], 'DisP_island'] = 'island'
    df_island[["chr", "start", "end", "signal", "rank", "diff", "DisP_island"]].to_csv(options['--output'] + '.island.txt', sep="\t", index = False)

    if options['--plot']:
        print('Hockey plot island signals vs rank...')
        sns.relplot(data = df_island, x = "rank", y = "signal", hue = "DisP_island", s = 10, edgecolor = "none")
        plt.xlabel("DisP-seq merged signals rank")
        plt.ylabel("DisP-seq signal")
        plt.title("DisP islands")
        plt.savefig(options['--output'] + '.pdf')
    file2rm = '%s,%s' % (options['--output'] + '.merged.peak', options['--output'] + '.sort.peak')
    for f in file2rm.strip().split(","):
            os.remove(f)
