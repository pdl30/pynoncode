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
from pynoncode import run_bowtie
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

def combine_samfiles(multi=False, clipped=False):
	#Seperate out clipped and unclipped!
	#Look at naming!
	if multi:
		sam1 = "unclipped_multimap.sam"
		sam2 = "clipped_multimap.sam"
		bam1 = "unclipped_multimap.bam"
		bam2 = "clipped_multimap.bam"
		out = open("multi_mapped.sam", "w")
	else:
		sam1 = "unclipped_unique.sam"
		sam2 = "clipped_unique.sam"
		bam1 = "unclipped_unique.bam"
		bam2 = "clipped_unique.bam"
		out = open("unique_mapped.sam", "w")
	#Convert unclipped sam to bam

	#Converts sam to bam
	bam1_o = open(bam1, "w")
	a = pysam.view("-bS", sam1)
	for r in a:                                     
		bam1_o.write(r)
	bam1_o.close()
	#Converts clipped sam to bam
	if clipped == True:
		if os.stat(sam2).st_size > 0: #Checking file is not empty
			try:
				bam2_o = open(bam2, "w")
				b = pysam.view("-bS", sam2)
				for r in b:                                     
					bam2_o.write(r)
				bam2_o.close()
			except:
				print "Samtools raised error, will assume Sam file is empty!"
			#Merge clipped and unclipped
			input_filenames = ["-f", bam1, bam2]
			output_filename = "tmp1.bam"
			merge_parameters = [output_filename] + input_filenames
			pysam.merge(*merge_parameters)
			pysam.sort("-n", "tmp1.bam", "tmp2" )
			subprocess.call(["rm", sam2, bam2])
	else:
		#If no clipped bam, just sort 
		pysam.sort("-n", bam1, "tmp2" )
	#Converts file to sam
	d = pysam.view("-h", "tmp2.bam")
	for r in d:                                     
		out.write(r)
	subprocess.call(["rm", "tmp2.bam", "tmp1.bam", sam1, bam1])

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

def combine_reports():
	txtlist = [ f for f in os.listdir(".") if f.endswith("_report.txt") ]
	output = open('mapping_report.txt', "w")
	for txt in txtlist:
		output.write("{}\n".format(txt))
		with open(txt) as f:
			for line in f:
				line = line.rstrip()
				if line.startswith("Warning"):
					pass
				else:
					output.write("{}\n".format(line)),
		os.remove(txt)
	output.close()

def main():	
	parser = argparse.ArgumentParser(description='Processes ncRNA samples from fastq files to sam file.\n Please ensure FASTQ files are in current directory.\n ')
	parser.add_argument('-f', '--fastq', help='Single end fastq', required=False)
	parser.add_argument('-p', '--paired', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	parser.add_argument('-i', '--index', help='Path to bowtie1 index', required=True)
	parser.add_argument('-c', action='store_true', help='Clip tRNA ends', required=False)
	parser.add_argument('-o', '--outdir', help='Output results to this directory', required=True)
	args = vars(parser.parse_args())
	index = args["index"]
	outdir = args["outdir"]
	path = os.getcwd()
	outpath = os.path.join(args["outdir"], "pynoncode")
	if os.path.isdir(outdir):
		print "Output directory exists"
	else:
		subprocess.call(["mkdir", outdir])
		subprocess.call(["mkdir", outpath])

	if args["paired"]:
		fq1 = args["paired"][0]
		fq2 = args["paired"][1]

		print("\nCreating Fasta files...\n"),
		parse_paired_fastq(fq1, fq2, outpath) 
		os.chdir(outpath)

		print("\nRunning Bowtie...\n"),
		run_bowtie.paired_bowtie(index)
		
		if args["c"]:
			strip_ends(True)
			run_bowtie.paired_bowtie(index, True)
			combine_samfiles(clipped=True) 
			combine_samfiles(multi=True, clipped=True) #Multimapper combination
		else:
			combine_samfiles() #Can't parallelise from here on!
			combine_samfiles(multi=True) #Multimapper combination

	elif args["fastq"]:
		fq1 = args["fastq"]

		print("\nCreating Fasta files...\n"),
		parse_single_fastq(fq1, outpath)
		os.chdir(outpath)

		print("\nRunning Bowtie...\n"),
		run_bowtie.single_bowtie(index)
		if args["c"]:
			strip_ends(False)
			run_bowtie.single_bowtie(index, True)
			combine_samfiles(clipped=True) #Can't parallelise from here on!
			combine_samfiles(multi=True, clipped=True) #Multimapper combination
		else:
			combine_samfiles()
			combine_samfiles(multi=True) 
	combine_reports()
	falist = [ f for f in os.listdir(".") if f.endswith(".fa") ]
	for f in falist:
		os.remove(f)
