#!/usr/bin/python

########################################################################
# 20 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import argparse
import subprocess
import sys, re, os
import pybedtools
import pkg_resources

def create_full_bedfile(ens):
	output = open("full_bed_file.BED", "w")
	with open("fragment_counts.txt") as f:
		for line in f:
			line = line.rstrip()
			word=  line.split("\t")
			score = int(float(word[4]))
			if ens:
				m = re.match("MT", word[0])
				if m:
					chrom = "chrM"
				else:
					chrom = "chr"+ word[0]
			else:
				chrom = word[0]
			for i in range(score):
				output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, word[1], word[2], word[3], word[4], word[5])),
	output.close()

def convertBed_bigWig(genome, chromsizes):
	inbed = pybedtools.BedTool("full_bed_file.BED")
	sortbed = inbed.sort()
	outcov = sortbed.genome_coverage(bg=True, genome=genome)
	outcov.saveas("pynoncode.bedGraph")
	command = ["bedGraphToBigWig", "pynoncode.bedGraph", chromsizes, "pynoncode.bw"]
	subprocess.call(command)

def cleanup():
	os.remove("pynoncode.bedGraph")
	os.remove("full_bed_file.BED")

def main():
	parser = argparse.ArgumentParser(description='Processes an ncpipe processed folder to create bigwigs\n')
	parser.add_argument('-i','--input', help='Ncpipe formatted directory', required=True)
	parser.add_argument('-p', help='Are samples paired end?', action='store_true', required=False)
	parser.add_argument('-e', help='Is sample aligned to Ensembl formatted genome?', action='store_true', required=False)
	parser.add_argument('-g', '--genome', help='Options are mm10/hg19', required=True)
	args = vars(parser.parse_args())
	chrom = pkg_resources.resource_filename('pynoncode', 'data/{}.chrom.sizes'.format(args["genome"]))
	os.chdir(args["input"])
	create_full_bedfile(args["e"])
	convertBed_bigWig(args["genome"], chrom)
	cleanup()
