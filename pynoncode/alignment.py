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
from pynoncode import run_bowtie, fasta_parsers
from itertools import izip

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

	samfile = pysam.Samfile(sam, "r")
	if paired:	
		first_reads = defaultdict(list)
		second_reads = defaultdict(list)
		for read in samfile.fetch() :
			if read.tid == -1:
				pass
			else:
				if read.is_proper_pair and read.is_read1:
					first_reads[read.qname].append(read)
				elif read.is_proper_pair and read.is_read2:
					second_reads[read.qname].append(read)
		for ids in first_reads:
			for f_reads in first_reads[ids]: #Loops over all those with the same IDs. Shouldn't matter if they are multiple mapped
				aligned = m.get(((f_reads.qname, samfile.getrname(f_reads.tid), f_reads.pos+1)), None)
				if aligned != "no_feature": #Exclude those uninteresting ones

					s_reads = second_reads[ids] #List of second reads
					for s_read in s_reads:
						if f_reads.pnext == s_read.pos: #If position of next read is same as next reads start positions
							paired_read = s_read
					strand = "+"
					if f_reads.is_reverse:
						strand = "-"
					strand2 = "+"
					if s_read.is_reverse:
						strand2 = "-"
					bedout.write("{}\t{}\t{}\t{}\t{}\t{}\t1\t{}\n".format(samfile.getrname(f_reads.tid), f_reads.pos, f_reads.aend, f_reads.qname, f_reads.seq, strand, aligned)),
					bedout.write("{}\t{}\t{}\t{}\t{}\t{}\t2\t{}\n".format(samfile.getrname(paired_read.tid), paired_read.pos, paired_read.aend, paired_read.qname, paired_read.seq, strand2, aligned)),
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

def convert_and_sort(sam):
	#No need to do this now!
	name = re.sub(".sam$", "", sam)
	command1 = "samtools view -bS {0}.sam > {0}.bam\n".format(name)
	command2 = "samtools sort -n {0}.bam {0}_sort\n".format(name)
	command3 = "samtools view -h -o {0}_sort.sam {0}_sort.bam\n".format(name)
	subprocess.call(command1, shell=True)
	subprocess.call(command2, shell=True)
	subprocess.call(command3, shell=True)
	os.remove("{0}.bam".format(name))
	os.remove("{0}_sort.bam".format(name))

def main():	
	parser = argparse.ArgumentParser(description='Processes ncRNA samples from fastq files to sam file.\n Please ensure FASTQ files are in current directory.\n ')
	parser.add_argument('-f', '--fastq', help='Single end fastq', required=False)
	parser.add_argument('-p', '--paired', help='Paired end fastqs. Please put them in order!', required=False, nargs='+')
	parser.add_argument('-i', '--index', help='Path to bowtie1 index', required=True)
	parser.add_argument('-c', action='store_true', help='Clip tRNA ends', required=False)
	parser.add_argument('-o', '--outdir', help='Output results to this directory', required=True)
	parser.add_argument('-a', '--gtf', help='Path to ncRNA GTF file. If not provided, will use packages mouse ensembl formatted GTF file', required=False)
	parser.add_argument('-g', '--genome', help='Genome samples are aligned to. Options are hg19/mm10', required=True)
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

	if args["gtf"]:
		gtf = args["gtf"]
	else:
		gtf = pkg_resources.resource_filename('pynoncode', 'data/{}_ncRNA.gtf'.format(args["genome"]))

	if args["paired"]:
		fq1 = args["paired"][0]
		fq2 = args["paired"][1]

		print("\nCreating Fasta files...\n"),
		fasta_parsers.parse_paired_fastq(fq1, fq2, outpath) 
		os.chdir(outpath)

		print("\nRunning Bowtie...\n"),
		run_bowtie.paired_bowtie(index)
		
		if args["c"]:
			fasta_parsers.strip_ends(True)
			run_bowtie.paired_bowtie(index, True)
			combine_samfiles(clipped=True) 
			combine_samfiles(multi=True, clipped=True) #Multimapper combination
		else:
			combine_samfiles() #Can't parallelise from here on!
			combine_samfiles(multi=True) #Multimapper combination

	elif args["fastq"]:
		fq1 = args["fastq"]

		print("\nCreating Fasta files...\n"),
		fasta_parsers.parse_single_fastq(fq1, outpath)
		os.chdir(outpath)

		print("\nRunning Bowtie...\n"),
		run_bowtie.single_bowtie(index)
		if args["c"]:
			fasta_parsers.strip_ends(False)
			run_bowtie.single_bowtie(index, True)
			combine_samfiles(clipped=True) #Can't parallelise from here on!
			combine_samfiles(multi=True, clipped=True) #Multimapper combination
		else:
			combine_samfiles()
			combine_samfiles(multi=True) 
	combine_reports()
	
	#Annotation part of the script
	annotate_sam("unique_mapped.sam", gtf)
	annotate_sam("multi_mapped.sam", gtf)
	convert_and_sort("unique_mapped.sam")
	convert_and_sort("multi_mapped.sam")
	convert_sam_bed("unique_mapped_sort.sam", "unique_mapped.samout", args["paired"], "unique_mapped.BED")
	convert_sam_bed("multi_mapped_sort.sam", "multi_mapped.samout", args["paired"], "multi_mapped.BED")
	cleanup()
	useful_fastas = ["original_fasta_1.fa", "original_fasta_2.fa", "original_fasta.fa", "clipped_1.fa", "clipped_2.fa", "clipped_fasta.fa"]
	falist = [ f for f in os.listdir(".") if f.endswith(".fa") ]
	for f in falist:
		if f not in useful_fastas:
			os.remove(f)
