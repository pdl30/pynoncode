#!/usr/bin/python

########################################################################
# 28 Apr 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import argparse
import subprocess
import sys, re, os
from collections import defaultdict
import pysam
import pybedtools
from logging import info
import ncpipe
from itertools import izip
import pkg_resources


def convert_and_sort(sam):
	name = re.sub(".sam$", "", sam)
	command1 = "samtools view -bS {0}.sam > {0}.bam\n".format(name)
	command2 = "samtools sort {0}.bam {0}_sort\n".format(name)
	command3 = "samtools index {0}_sort.bam\n".format(name)
	subprocess.call(command1, shell=True)
	subprocess.call(command2, shell=True)
	subprocess.call(command3, shell=True)
#	os.remove(sam)
	os.remove("{0}.bam".format(name))

def convert_sam_bed(sam, samout, paired, bed):
	m = {}
	with open(samout) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if len(word) > 13:
				mapped = word[14].split(":")
				mapped_to = re.sub("__", "", mapped[2])
				#Read name, chromosome, start - 1 based
				m[word[0], word[2], int(word[3])] = mapped_to
	bedout = open(bed, "w")

	samfile = pysam.Samfile(sam, "rb")
	if paired:	
		for read in samfile.fetch() :
			if read.tid == -1:
				pass
			else:
				if read.is_proper_pair and read.is_read1 and read.tid == read.rnext:
					pointer = samfile.tell() # pointer to the current position in the BAM file
					try: 
						mate = samfile.mate(read) 
					except ValueError: 
						continue 
					finally: 
						 samfile.seek(pointer) 
					strand = '+'
					if read.is_reverse:
						strand = '-'
					aligned = m.get(((read.qname, samfile.getrname(read.tid), read.pos+1)), None)
					bedout.write("{}\t{}\t{}\t{}\t{}\t{}\t1\t{}\n".format(samfile.getrname(read.tid), read.pos, read.aend, read.qname, read.seq, strand, aligned)),
				#elif read.is_proper_pair and read.is_read2:
					strand2 = '+'
					if mate.is_reverse:
						strand2 = '-'
					aligned = m.get(((mate.qname, samfile.getrname(mate.tid), mate.pos+1)), None)
					bedout.write("{}\t{}\t{}\t{}\t{}\t{}\t2\t{}\n".format(samfile.getrname(mate.tid), mate.pos, mate.aend, mate.qname, mate.seq, strand2, aligned)),
	else:
		for read in samfile.fetch() :
			if read.tid == -1:
				pass
			else:
				strand = '+'
				if read.is_reverse:
					strand = '-'
				aligned = m.get(((read.qname, samfile.getrname(read.tid), read.pos+1)), None)
				bedout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(samfile.getrname(read.tid), read.pos, read.aend, read.qname, read.seq, strand, aligned)),
	bedout.close()
	samfile.close()
	

def annotate_sam(sam_file, gtf_file, nctype=None):
	name = re.sub(".sam", "", sam_file)
	if nctype == None:
		command = ["htseq-count", "--quiet", "--samout={}.samout".format(name), "--minaqual=0", "--idattr=transcript_id", "{}.sam".format(name),  gtf_file]
	else:
		if nctype == "miRNA" or nctype == "tRNA":
			command = ["htseq-count", "--quiet", "--type={}".format(nctype), "--minaqual=0", "--samout={}.samout".format(name), "{}.sam".format(name), gtf_file]
		else:
			raise Exception("Unrecognised ncRNA type!")
	with open(os.devnull, 'w') as devnull:
		subprocess.call(command, stdout=devnull)

def cleanup():
	os.remove("unique_mapped.samout")
	os.remove("multi_mapped.samout")

def main():
	parser = argparse.ArgumentParser(description='Processes ncRNA samples from fastq files to sam file.\n Please ensure FASTQ files are in current directory.\n ')
	parser.add_argument('-i', '--input', help='Input directory after ncalign has been run', required=True)
	parser.add_argument('-g', '--gtf', help='Path to ncRNA GTF file. If not provided, will use packages mouse ensembl formatted GTF file', required=False)
	#parser.add_argument('-type', help='Type of ncRNA to count, if not used, all ncRNAs are considered', required=False) #not implemented yet
	parser.add_argument('-p', '--paired', help='Experiment is paired end', action="store_true", required=False)
	args = vars(parser.parse_args())
	if args["gtf"]:
		gtf = args["gtf"]
	else:
		gtf = pkg_resources.resource_filename('pynoncode', 'data/mm10_ncRNA.gtf')

	path = os.path.join(args["input"], "pynoncode")
	os.chdir(path)
	print("Annotating Sam File...\n"),
	annotate_sam("unique_mapped.sam", gtf)
	annotate_sam("multi_mapped.sam", gtf)
	convert_and_sort("unique_mapped.sam")
	convert_and_sort("multi_mapped.sam")
	convert_sam_bed("unique_mapped_sort.bam", "unique_mapped.samout", args["paired"], "unique_mapped.BED")
	convert_sam_bed("multi_mapped_sort.bam", "multi_mapped.samout", args["paired"], "multi_mapped.BED")
	cleanup()