#pyrnatools 

### Installation

Clone this repository and then:

```bash
$ cd pynoncode/
$ python setup.py install --user
```

This will install the scripts found in the pyrnatools/scripts directory. For more information on the individual scripts, use the --help command after each script. 

##Dependencies

- bowtie version 1.0.0 or greater
- HTseq-count version 0.6.0 or greater
- samtools Version: 0.1.19 or greater
- bedGraphToBigWig from UCSC binaries

##Core Pipeline
- pynon_align.py - Converts FASTQ to FASTA and runs bowtie aligner to exact both unique and multimapped sequences
- pynon_anno.py - Annotates SAM files and converts them to BED format
- pynon_count.py - Creates transcript counts and fragment counts files from BED files using if necessary multiple mapped reads
- pynon_viz.py - Converts fragments table to UCSC formatted bigWigs
