#pynoncode 

### Installation

Clone this repository and then:

```bash
$ cd pynoncode/
$ python setup.py install
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
- bowtie version 1. [Link](http://bowtie-bio.sourceforge.net/index.shtml)
- samtools version: 0.1.19 or greater. [Link](http://www.htslib.org/) 
- A [R](http://www.r-project.org/) installation.
- DESeq2 [Link](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html).
- The remaining python dependencies should be installed automatically. These are:
  - pybedtools [Link](http://pythonhosted.org/pybedtools/)
  - pysam [Link](https://github.com/pysam-developers/pysam)
  - HTSeq [Link](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)

##Core Pipeline

- pynon_align.py:
  Please trim adapters before using this program. Some programs for doing this are [cutadapt](https://code.google.com/p/cutadapt/). 
  This converts FASTQ to FASTA and runs bowtie aligner to exact both unique and multimapped sequences. If the -c option is specified, pynon_align.py will trim tRNA CCA ends from the unmapped reads and then rerun the alignment step. The program will then use a GTF file to annotate the mapped fragments using HTSeq. The resulting files can be found in pynoncode directory and are called unique_mapped.BED and multi_mapped.BED

- pynon_count.py:
  This uses the uniquely mapped reads and if -m is specified, will use the multi-mapped reads to create counts files of both the transcripts and fragments. 
  Multi-mapped reads are distributed according to the fraction of unique counts the transcripts have.

- pynon_ucsc.py:
  This create a UCSC formatted bigWig for use on the UCSC genome browser. The bigWig is called pynoncode.bw and can be found in the results directory.

- pynon_diff.py:
  This uses DESeq2 to perform differential expression on transcripts and fragments. 
  For examples of configuration files please see [here](https://github.com/pdl30/pynoncode/tree/master/configuration_examples)

- pynon_report.py:
  Plots profiles of transcripts and creates a HTML report of the differentially expressed fragments and transcripts. 
  For examples of configuration files please see [here](https://github.com/pdl30/pynoncode/tree/master/configuration_examples)

##Annotation

- Inlcuded in this package are a mouse and human noncoding GTF. For more information, see [here](https://github.com/pdl30/pynoncode/tree/master/pynoncode/data)
- However if you wish to use your own annotation, please make sure it is in GTF format. For more information see [here](http://www.ensembl.org/info/website/upload/gff.html).
- Also please note the chromosome names in your GTF file and use the options to convert it to UCSC notation if required. 
