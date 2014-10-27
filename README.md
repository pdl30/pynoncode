#ncRNA Analysis Toolkit

To use the package, please download a copy, unzip and then add the path to your PATH. For example: 

export PATH=$PATH:~/Programs/ncrna_toolkit/

Also add the path to the downloaded directory to your ~/.bashrc for the python modules. For example:

export PYTHONPATH=$PYTHONPATH:~/Programs/

These Python packages will be installed on first usage:

pysam

pybedtools

sqlite3

###R Packages required. These need to be manually installed:

DESeq

RSQLite


#Usage:

##The pipeline consists of 3 scripts:



###fragment_counts.py -h

usage: fragment_count.py [-h] [-f FASTQ] [-p PAIRED [PAIRED ...]] -i INDEX
                         [-g GTF] -o OUTDIR [-t TYPE]

Processes ncRNA samples from fastq files to sam file. Please ensure FASTQ files are in current directory.

Arguments:

  -h, --help            show this help message and exit

  -f FASTQ, --FASTQ FASTQ Single end fastq

  -p PAIRED [PAIRED ...], --PAIRED PAIRED [PAIRED ...] Paired end fastqs. Please put them in order!

  -i INDEX, --INDEX INDEX Path to bowtie2 index

  -g GTF, --GTF GTF     Path to ncRNA GTF file. If not provided, will use packages GTF file

  -o OUTDIR, --OUTDIR OUTDIR Output results to this directory

  -t TYPE, --TYPE TYPE  Type of ncRNA to count, if not used, all ncRNAs are considered


The options -p and -f are mutually exclusive depending on the type of Fastq's.    





###diff_fragments.py

usage: diff_fragments.py [-h] -d DESIGN [-p] [-o OUTDIR]

Processes ncRNA counts.

Arguments:

  -h, --help            show this help message and exit

  -d DESIGN, --DESIGN DESIGN Text file containing sample names and conditions

  -p, --PAIRED          Flag indicating samples are paired end

  -o OUTDIR, --OUTDIR OUTDIR Output results to this directory, don't use! not implemented yet


This creates a Fragments sqlite database in current directory which will be quiered for the Fragments data. 

Do not run multiple instances of this script in the same directory!

A design matrix is required for usage. This must be a tab-deliminated file containing a header row! 

Only 2 conditions are supported at this time! The first column must correspond to a directory containing the files outputted by fragment_counts.py!

An example is:


Sample Condition

do2703  WT

do2704  WT

do2707  KO

do2708  KO

This script will create a file called run_diff_rcode.R in the current directory. If the differential analysis fails, the R code used will be present in this file for debugging and rerunning!

The most common error at this stage is:

Error in parametricDispersionFit(means, disps) :

If this occurs, insert fitType = "local" at estimateDispersions line and call from the command line:

Rscript run_diff_rcode.R



###samtobigWig.py

usage: samtobigWig.py [-h] -i DIR -p -c CHROM [-o OUT] [-s SPLIT]

Processes ncRNA sam files to bigWig files. Please ensure Sam files are in
results directory

Arguments:

  -h, --help            show this help message and exit

  -i DIR, --DIR DIR     The path to the directory processed by fragment_count.py. Results will be outputted in this directory

  -p, --PAIRED          Flag indicating paired end samples

  -c CHROM, --CHROM CHROM Chromosome sizes. This must correspond to the index used by bowtie2!

  -o OUT, --OUT OUT     Name of final bigwig files, not implemented yet!

  -s SPLIT, --SPLIT SPLIT Split the resulting bigwig files by strand?, not implemented yet!


Creates BigWigs in specified directory! 
Warning! Creates very large files in specified directory can be deleted after usage.