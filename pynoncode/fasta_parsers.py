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
from itertools import izip

def parse_paired_fastq(fq1, fq2, outdir):
	dict2 = defaultdict(int)
	count_dict = defaultdict(int)
	f1=open(fq1)
	f2=open(fq2)
	for line1, line2 in izip(f1, f2):
		line1 = line1.rstrip()
		id1 = line1.split("#")
		line2 = line2.rstrip()
		id2 = line2.split("#")
		try:
			id11 = next(f1)
			read1 = id11.rstrip()
			id22 = next(f2)
			read2 = id22.rstrip()
			reads = "{}\t{}".format(read1, read2)
			dict2[reads] += 1
			crap = next(f1)
			crap2 = next(f1)
			crap = next(f2)
			crap2 = next(f2)
		except StopIteration:
			break
	seen = {}
	name1 = "original_fasta_1.fa"
	name2 = "original_fasta_2.fa"
	count = 1
	output1 = open(outdir + "/" + name1, "wb")
	output2 = open(outdir + "/" + name2, "wb")
	for key in dict2.keys():
		reads = key.split("\t")
		output1.write(">ID:{}\n{}\n".format(count, reads[0])),
		output2.write(">ID:{}\n{}\n".format(count, reads[1])),
		count_dict[count] = dict2[key]
		count += 1
	output3 = open(outdir + "/" + "count_dict.txt", "w")
	for key in count_dict.keys():
		output3.write("{}\t{}\n".format(key, count_dict[key])),

def parse_single_fastq(fq1, outdir):
	dict2 = defaultdict(int)
	count_dict = defaultdict(int)
	f1=open(fq1)
	for line1 in f1:
		line1 = line1.rstrip()
		id1 = line1.split("#")
		try:
			id11 = next(f1)
			read1 = id11.rstrip()
			dict2[read1] += 1
			crap = next(f1)
			crap2 = next(f1)
		except StopIteration:
			break
	seen = {}
	name1 = "original_fasta.fa"
	count = 1
	output1 = open(outdir + "/" + name1, "wb")
	for key in dict2.keys():
		reads = key.split("\t")
		output1.write(">ID:{}\n{}\n".format(count, reads[0])),
		count_dict[count] = dict2[key]
		count += 1
	output3 = open(outdir + "/" + "count_dict.txt", "w")
	for key in count_dict.keys():
		output3.write("{}\t{}\n".format(key, count_dict[key])),

def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

	
def stripper(fasta):
	result = {}
	with open(fasta) as f:
		for name, seq in read_fasta(f):
			bases = list(seq)
			end1 = bases[-3:]
			end1 = ''.join(end1)
			if end1 == "CCA":
				tmpseq = bases[:-3]
				seq = ''.join(tmpseq)
			end2 = bases[-4:]
			end2 = ''.join(end2)
			if end2 == "CCAC":
				tmpseq = bases[:-4]
				seq = ''.join(tmpseq)
			end3 = bases[-5:]
			end3 = ''.join(end3)
			if end3 == "CCACC":
				tmpseq = bases[:-5]
				seq = ''.join(tmpseq)
			end4 = bases[-6:]
			end4 = ''.join(end4)
			if end4 == "CCACCA":
				tmpseq = bases[:-6]
				seq = ''.join(tmpseq)
			result[name] = seq
	return result

def strip_ends(paired):
	if paired == True:
		output1 = open("clipped_1.fa", "w")
		output2 = open("clipped_2.fa", "w")
		
		data1 = stripper("unclipped_multi_unmapped_1.fa")
		data2 = stripper("unclipped_multi_unmapped_2.fa")
		for key in sorted(data1.keys()):
			output1.write("{}\n{}\n".format(key, data1[key])),
		for key in sorted(data2.keys()):
			output2.write("{}\n{}\n".format(key, data2[key])),
	else:
		data1 = stripper("unclipped_multi_unmapped.fa")
		output1 = open("clipped_fasta.fa", "w")

		for key in sorted(data1.keys()):
			output1.write("{}\n{}\n".format(key, data1[key])),