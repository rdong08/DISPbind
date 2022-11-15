# *align*

`DISPbind align` aligns DisP-seq with BWA.

## Usage and option summary

### Usage:

```
 DISPbind.py align [options] -i INDEX -a FQ1 -b FQ2 -o OUT -n NAME
```

### Options:

```
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

## Notes about options

1. The -o OUT indicate the output folder and -n NAME is the prefix of output bam files.
2. It will overwrite the output directory, so please double check when setting the utput directory.

## Output

`DISPbind align` will create a folder with bam file(`name.sorted.deduped.bam`) and bigwig file(`name.bw`). The Bam would be used for DisP-seq peak calling and bigwig file will be used for further DisP island annotation.

```
output/name.sorted.deduped.bam
output/name.bw
```
