# DISPbind: Disorder protein genomic binding analysis toolkit for DisP-seq

A schematic overview of DisP-seq
-----------------------------------
<img src="https://github.com/rdong08/DISPbind/blob/main/docs/image/DisP-seq.png" width="1000">

## Features

* Genome mapping and bigwig generation of DisP-seq sequences
* Identification of DisP island 

## Prerequisites
* numpy (https://numpy.org/)
* pysam (https://pysam.readthedocs.io/en/latest/api.html
* pybedtools (https://daler.github.io/pybedtools/)
* pandas (https://pandas.pydata.org/)
* pyBigWig (https://github.com/deeptools/pyBigWig)
* docopt (http://docopt.org/)
* matplotlib (https://matplotlib.org/)
* seaborn (https://seaborn.pydata.org/)
* MACS (https://github.com/macs3-project/MACS)

## Installation

### Install by using pip
```bash
pip install DISPbind
```

## Usage and example
-----
### Step1: Genome mapping and bigwig generation of DisP-seq sequences
<img src="https://github.com/rdong08/DISPbind/blob/main/docs/image/mapping.png" width="1000">

```
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
```

```bash
DISPbind align -i bwa_index_hg19/hg19.fa -n test -a test_file/test_SKNMC_bisox_rep_R1.fastq -b test_file/test_SKNMC_bisox_rep_R2.fastq -o test_out -p 1 -g hg19.chrom.sizes
```

Transfer Bam file to bigwig file **(Only if you do not run the align step and could provide DisP-seq Bam file)**
```
Usage: DISPbind bam2bw [options] -b bam -n NAME -o OUT

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -m MQ --mquality=MQUALITY      Mapping quality. [default: 10]
    -g GSIZE --gsize=GSIZE         Genome size file.
    -n NAME --name=NAME            Output file name. [default: bwa_out]
    -b BAM --bam=BAM               Input BAM file.
    -o OUT --output=OUT            Output directory. [default: alignment]
```

```bash
DISPbind bam2bw -b test/test.bam -n test -o test_out -g hg19.chrom.sizes
```

### Step2: DisP-seq peak calling by MACS
```bash
macs2 callpeak --nomodel -B --SPMR -f BAMPE -g hs -t Bisox.sorted.deduped.bam -c Input.sorted.deduped.bam --outdir callpeak/ -n bisox_peaks --broad
```

### Step3: Identify DisP islands from bigwig and peak file
<img src="https://github.com/rdong08/DISPbind/blob/main/docs/image/island.png" width="1000">

```
Usage: DISPbind island [options] -p PEAK (-b BW | -l LIST) -o OUT

Options:
    -h --help                             Show help message.
    -v --version                          Show version.
    -p PEAK --peak=PEAK                   DisP peaks.
    -b BW --bigwig=BW                     Bigwig file for corresponding sample.
    -o OUTPREFIX --output=OUTPREFIX       Output island file. [default: DisP]
    -l LIST --list=LIST                   Bigwig replicate (will calculate average signals for replicates)
    --plot                                Hockey plot for signal vs rank.
```

```bash
DISPbind island -o island.test -p bisox_peaks.bed -b bisox.sorted.fragments.bw --plot
```

Format of output `island.test.island`:

| Field       | Description                        |
| :---------: | :----------------------------------|
| Chrom       | Chromosome                         |
| Start       | Start of DisP merged region        |
| End         | End of DisP merged region          |
| Signal      | Sum of DisP signal in merged region|
| Rank        | Rank of DisP signal                |
| Island      | Label Island or not Island         |


## Citation

Yu-Hang Xing\*, Rui Dong\*, Lukuo Lee, Shruthi Rengarajan, Nicolò Riggi, Gaylor Boulay and Miguel N. Rivera#. **DisP-seq reveals the genome-wide functional organization of DNA-associated disordered proteins** *Under review*

## License
Copyright (C) 2022 Rivera Lab. See [LICENSE](about/license.md) for license rights and limitations (MIT).
