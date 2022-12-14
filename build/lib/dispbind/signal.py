import pyBigWig
import numpy as np
import argparse
import re
import os

def fetch_genome_size(genome_size_file_in):
    gs = dict()
    for line in open(genome_size_file_in):
        cols = line.strip().split()
        gs[cols[0]] = int(cols[1])
    return(gs)

def fetch_bed_loc(bed_file_in):
    bed_loc_out = []
    for line in open(bed_file_in):
        cols = line.strip().split()
        loc_start = int(cols[1])
        loc_end = int(cols[2])
        if (re.match("chr",cols[0])):
            loc_chr = cols[0]
        else:
            loc_chr = str("chr")+str(cols[0])
        bed_loc_out.append(cols[0:3])
    return(bed_loc_out)

def fetch_signal(bw_file_in, bed_loc_in):
    bw_signal = []
    bw_parse = pyBigWig.open(bw_file_in)
    for cols in bed_loc_in:
        loc_start = int(cols[1])
        loc_end = int(cols[2])
        chr_info = cols[0]
        if (sorted(list(bw_parse.chroms().keys()))[0] == "1"):
            chr_info = re.sub(r'chr', '', cols[0])
        if (chr_info in bw_parse.chroms().keys()):
            vals = bw_parse.values(chr_info, loc_start, loc_end)
            vals = np.array(vals)
            vals[np.isnan(vals)] = 0
            bw_signal.append(np.sum(vals))
        else:
            bw_signal.append(0)
    bw_parse.close()
    return(bw_signal)
