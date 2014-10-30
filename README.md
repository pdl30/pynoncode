#pynoncode 

### Installation

Clone this repository and then:

```bash
$ cd pynoncode/
$ python setup.py install --user
```

This will install the scripts found in the pyrnatools/scripts directory. For more information on the individual scripts, use the --help command after each script. 

##Dependencies
- numpy - This can be installed using:
```bash
	pip install numpy
```
- cython - This can be installed using:
```bash
	pip install cython
```
- bowtie version 1.0.0 or greater. [Link](http://bowtie-bio.sourceforge.net/index.shtml)
- samtools version: 0.1.19 or greater. [Link](http://www.htslib.org/) 

##Core Pipeline
- pynon_align.py - Converts FASTQ to FASTA and runs bowtie aligner to exact both unique and multimapped sequences
- pynon_anno.py - Annotates SAM files and converts them to BED format
- pynon_count.py - Creates transcript counts and fragment counts files from BED files using if specified multiple mapped reads
- pynon_viz.py - Converts fragments table to UCSC formatted bigWigs

##Annotation
- Inlcuded in this package are a mouse and human GTF. For more information, see [here](https://github.com/pdl30/pynoncode/tree/master/pynoncode/data)
- However if you wish to use your own annotation, please make sure it is in GTF format. For more information see [here](http://www.ensembl.org/info/website/upload/gff.html).
- Also please note the chromosome names in your GTF file and use the options to convert it to UCSC notation if required. 
