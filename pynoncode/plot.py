#!/usr/bin/python

########################################################################
# 31 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, re, os
import ConfigParser
import itertools
import argparse
from collections import defaultdict
import pkg_resources
import HTSeq
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool, Manager
from pynoncode import web_templates

def read_custom_input(ifile):
	trans = {}
	with open(ifile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			trans[word[0]] = 1
	return trans

def read_trans_input(ifile, pval):
	trans = {}
	with open(ifile) as f:
		header = next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[5] == "NA":
				pass
			elif float(word[5]) <= pval:
				trans[word[0]] = (word[2], word[5], word[6]) #LFC Pvalue, Padj
	return trans

def read_frag_input(ifile, paired, pval):
	frags = {}
	with open(ifile) as f:
		header = next(f)
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if not paired:
				if word[5] == "NA":
					pass
				elif float(word[5]) <= float(pval):
					frags[word[0]] = (word[2], word[5], word[6])
			else:
				pairs = word[0].split("|")
				if word[5] == "NA":
					pass
				elif float(word[5]) <= pval:
					frags[pairs[0], pairs[1]] = (word[2], word[5], word[6]) #LFC. Pvalue, padj
	return frags ##Will start with single and then consider paired

def find_frag_transcripts(conditions, frags, paired):
	transcripts = {}
	for idir in conditions:
		with open(idir + "/fragment_counts.txt") as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if not paired: #This is adding every occurence, need to just find unique positions!
					if word[3] in frags:
						if word[6] not in transcripts:
							transcripts[word[6]] = {}
							transcripts[word[6]][word[0], word[1], word[2]] = (word[3], frags[word[3]])#Position of diff fragment. Make this list if more than one are involved???
						else:
							transcripts[word[6]][word[0], word[1], word[2]] = (word[3], frags[word[3]])
				else:
					if int(word[7]) == 1:
						next_line = next(f).rstrip()
						next_word = next_line.split("\t")
						if int(next_word[7]) == 2:
							if (word[3], next_word[3]) in frags:
								if word[6] not in transcripts:
									transcripts[word[6]] = {}
									transcripts[word[6]][word[0], word[1], word[2], next_word[1], word[2]] = (word[3], next_word[3], frags[word[3], next_word[3]]) #Contains both reads positions
								else:
									transcripts[word[6]][word[0], word[1], word[2], next_word[1], word[2]] = (word[3], next_word[3], frags[word[3], next_word[3]])
	return transcripts

#Could add region in plots!
def read_directories_for_transcripts(conditions, transcript_coords, paired):
	#Counts the incidence of transcripts in fragments file
	transcript_arrays = {} #Initialise this dictionary
	for trans in transcript_coords:
		transcript_arrays[trans] = {}
		for idir in conditions:
			length = int(transcript_coords[trans][2]) - int(transcript_coords[trans][1])
			transcript_arrays[trans][idir] = np.zeros(length, dtype='f')
	for idir in conditions:
		with open(idir + "/fragment_counts.txt") as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if not paired:
					if word[6] in transcript_coords:
						start = int(word[1]) - int(transcript_coords[word[6]][1])
						end = int(word[2]) - int(transcript_coords[word[6]][1])
						if start < 0:
							start = 0
						if end > int(transcript_coords[word[6]][2]):
							end = int(transcript_coords[word[6]][2])
						transcript_arrays[word[6]][idir][start:end] += float(word[4])#Add the count of that fragment to the transcript range it covers
				else:
					#I only need count once per pair! But I need to add this over all regions per pair!!
					if int(word[7]) == 1: #First pair
						if word[6] in transcript_coords: #Make sure transcript is important
							next_line = next(f).rstrip()
							next_word = next_line.split("\t")
							if int(next_word[7]) == 2: #Make sure they are properly paired
								p1_start = int(word[1]) - int(transcript_coords[word[6]][1])
								p1_end = int(word[2]) - int(transcript_coords[word[6]][1])
								if p1_start < 0:
									p1_start = 0
								if p1_end > int(transcript_coords[word[6]][2]):
									p1_end = int(transcript_coords[word[6]][2])
								transcript_arrays[word[6]][idir][p1_start:p1_end] += float(word[4]) #This is first pair

								p2_start = int(next_word[1]) - int(transcript_coords[word[6]][1]) #Second pair
								p2_end = int(next_word[2]) - int(transcript_coords[word[6]][1])
								if p2_start < 0:
									p2_start = 0
								if p2_end > int(transcript_coords[word[6]][2]):
									p2_end = int(transcript_coords[word[6]][2])
								transcript_arrays[word[6]][idir][p2_start:p2_end] += float(word[4]) #Must add same count to every pair

	return transcript_arrays

def invert_dict_nonunique(d):
    newdict = {}
    for k, v in d.iteritems():
        newdict.setdefault(v, []).append(k)
    return newdict

def average_arrays(conditions, transcript_arrays, transcript_coords):
	#Reverse the conditions dictionary and then average counts per conditions
	inv_conditions = invert_dict_nonunique(conditions)
	inv_array = {}
	for transcript in sorted(transcript_arrays):  #Naming is confusing, transcript here is an array
		length = int(transcript_coords[transcript][2]) - int(transcript_coords[transcript][1])
		inv_array[transcript] = {}
		for cond in inv_conditions:
			inv_array[transcript][cond] = np.zeros(length, dtype='f')
			count = 1
			for sample in inv_conditions[cond]:
				inv_array[transcript][cond] += transcript_arrays[transcript][sample]
				count += 1
			inv_array[transcript][cond] /= count
	return inv_array

def plot_trans_arrays(conditions, transcript_arrays, outputdir):
	#Plot sererately per transcript
	for transcript in sorted(transcript_arrays):
		c = 1
		for sample in sorted(transcript_arrays[transcript]):
			c +=1
			length_of_transcript = len(transcript_arrays[transcript][sample])
			base_label = np.array(xrange(length_of_transcript))
			plt.plot(base_label, transcript_arrays[transcript][sample], label="{}".format(sample)) 
		if c <= 4: #Control size of legend
			plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0., prop={'size':8})
		else:
			plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0., prop={'size':5})
			
		plt.ylabel('Read Count')
		plt.savefig(outputdir+'/plots/{}.png'.format(transcript))
		plt.close()

def plot_frag_arrays(conditions, transcript_arrays, outputdir, transcript_coords, transcripts_dict, paired):
	#Plot sererately per transcript
	for transcript in sorted(transcripts_dict): #Key is transcript, values are dict of positions and then fragments 
		a = 1
		for frag_pos in sorted(transcripts_dict[transcript]):
			c = 1
			for sample in sorted(transcript_arrays[transcript]):
				length_of_transcript = len(transcript_arrays[transcript][sample])
				base_label = np.array(xrange(length_of_transcript))

				c += 1 #Count number of samples for legend size	
				plt.plot(base_label, transcript_arrays[transcript][sample], label="{}".format(sample)) #Same as transcripts
			
			start_pos = int(frag_pos[1]) - int(transcript_coords[transcript][1])
			end_pos = int(frag_pos[2]) - int(transcript_coords[transcript][1])
			plt.axvspan(start_pos, end_pos, color='red', alpha=0.2)
			if paired:
				start_pos = int(frag_pos[3]) - int(transcript_coords[transcript][1])
				end_pos = int(frag_pos[4]) - int(transcript_coords[transcript][1])
				plt.axvspan(start_pos, end_pos, color='red', alpha=0.2)
		#Plot labels
			if c <= 4: #Control size of legend
				plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0., prop={'size':5})
			else:
				plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0., prop={'size':7})
			plt.ylabel('Read Count')
			plt.savefig(outputdir+'/plots/{}_{}.png'.format(transcript, a))
			plt.close()
			a += 1

#To reduce memory usage, just store interesting transcripts
def preprocess_gtf(gtf, transcripts):
	gtffile = HTSeq.GFF_Reader( gtf )
	exons = defaultdict(list)
	for feature in gtffile:
		if feature.type == "exon":
			if feature.attr["transcript_id"] in transcripts:
				exons[feature.attr["transcript_id"]].append(feature) #Just a list of exons for each transcript
	transcript_coords = {}
	for trans in sorted(exons):
		list_of_exons = exons[trans]
		if len(list_of_exons) == 1: #Don't care about exon numbers, just get start and end. Strand is unimportant
			transcript_coords[trans] = (exons[trans][0].iv.chrom, exons[trans][0].iv.start, exons[trans][0].iv.end, exons[trans][0].iv.strand) 
		else:
			#Strand is important, need to becareful here.
			exon_count = 1
			for exon in list_of_exons:
				if exon.iv.strand == "+":
					if exon.attr["exon_number"] == "1":
						chrom = exon.iv.chrom
						start = exon.iv.start
						strand = exon.iv.strand
					else:
						if exon.attr["exon_number"] > exon_count:
							end = exon.iv.end
				else: #Reverse start and end for negative strands
					if exon.attr["exon_number"] == "1":
						chrom = exon.iv.chrom
						end = exon.iv.end
						strand = exon.iv.strand
					else:
						if exon.attr["exon_number"] > exon_count:
							start = exon.iv.start
				exon_number = exon.attr["exon_number"]
			transcript_coords[trans] = (chrom, start, end, strand)
	return transcript_coords

def ConfigSectionMap(section, Config):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1


def write_trans_summary(transcripts, outdir):
	output = open(outdir + "/transcript_summary.tsv", "w")
	output.write("Transcript\tP-Value\tLFC\n"),
	for trans in sorted(transcripts):
		output.write("{}\t{}\t{}\n".format(trans, transcripts[trans][1], transcripts[trans][0])),
	output.close()
	
def write_frag_summary(transcripts, paired, outdir):
	output = open(outdir + "/fragment_summary.tsv", "w")
	if paired:
		output.write("Chromosome\tStart\tEnd\tRead1\tChromosome\tStart\tEnd\tRead2\tP Value\tLog Fold Change\tMapped Transcript\n"),
	else:
		output.write("Chromosome\tStart\tEnd\tRead1\tP Value\tLog Fold Change\tMapped Transcript\n"),
	for trans in sorted(transcripts):
		for frag_pos in sorted(transcripts[trans]):
			if paired:
				output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(frag_pos[0], frag_pos[1], frag_pos[2], transcripts[trans][frag_pos][0], frag_pos[0], frag_pos[3], frag_pos[4], 
					transcripts[trans][frag_pos][1], transcripts[trans][frag_pos][2][1], transcripts[trans][frag_pos][2][0], trans))
			else:
				output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(frag_pos[0], frag_pos[1], frag_pos[2], transcripts[trans][frag_pos][0], 
					transcripts[trans][frag_pos][1][1], transcripts[trans][frag_pos][1][0], trans))
	output.close()

def main():
	parser = argparse.ArgumentParser(description='Plots transcripts from pynon processed samples\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	#Parent Subparser
	base_subparser = argparse.ArgumentParser(add_help=False)
	base_subparser.add_argument('-c','--config', help='Config file, similar as the one supplied to pynon_diff.py. Please see documentation for more details', required=True)
	base_subparser.add_argument('-i','--input', help='Input file. Can be custom input file with transcript name on first column, transcripts or fragment files from DESEQ2', required=True)
	base_subparser.add_argument('-g','--gtf', help='GTF for annotation. If not supplied, will use the packages GTF')
	base_subparser.add_argument('-n','--genome', help='Sample genome name, options are hg19/mm10', required=True)
	base_subparser.add_argument('-p', help='Use if samples are paired end', action="store_true", required=False)
	base_subparser.add_argument('-a', help='Average sample according to sample conditions',  action="store_true", required=False)
	
	base_subparser.add_argument('-o','--outdir', help='Output directory', required=True)

	custom_parser = subparsers.add_parser('custom', help="Provide custom input file", parents=[base_subparser])
	trans_parser = subparsers.add_parser('trans', help="Provide differential transcript output from pynon_diff.py", parents=[base_subparser])
	trans_parser.add_argument('-v', '--pval', help='Pvalue cutoff for significance, default=0.1', default=0.1, required=False)
	frags_parser = subparsers.add_parser('frags', help="Provide differential fragment output from pynon_diff.py", parents=[base_subparser])
	frags_parser.add_argument('-v', '--pval', help='Pvalue cutoff for significance, default=0.1', default=0.1, required=False)
	args = vars(parser.parse_args())

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	conditions = ConfigSectionMap("Conditions", Config)

	if os.path.isdir(args["outdir"]):
		print "Output directory already exists, may overwrite existing files"
	else:
		os.mkdir(args["outdir"])
		os.mkdir("{}/plots".format(args["outdir"]))

	if args["gtf"]:
		gtf = args["gtf"]
	else:
		gtf = pkg_resources.resource_filename('pynoncode', 'data/{}_ncRNA.gtf'.format(args["genome"]))

	path_to_stuff = pkg_resources.resource_filename('pynoncode', 'data/') #Used for web templates css features

	if args["subparser_name"] == "custom": 
		transcripts = read_custom_input(args["input"])
		transcript_coords = preprocess_gtf(gtf, transcripts) #Annotation of transcripts
		transcript_arrays = read_directories_for_transcripts(conditions, transcript_coords, args["p"]) #Dict of numpy array containing counts of transcripts
		if args["a"]: #Average over conditions by reversing numpy array dict
			averaged_array = average_arrays(conditions, transcript_arrays, transcript_coords)
			plot_trans_arrays(conditions, averaged_array, args["outdir"])
		else:
			plot_trans_arrays(conditions, transcript_arrays, args["outdir"])

	elif args["subparser_name"] == "trans": 
		transcripts = read_trans_input(args["input"], args["pval"]) #Now contains LFC, Pvalue and padj
		transcript_coords = preprocess_gtf(gtf, transcripts) #Annotation of transcripts
		transcript_arrays = read_directories_for_transcripts(conditions, transcript_coords, args["p"]) #Dict of numpy array containing counts of transcripts
		if args["a"]: #Average over conditions by reversing numpy array dict
			averaged_array = average_arrays(conditions, transcript_arrays, transcript_coords)
			plot_trans_arrays(conditions, averaged_array, args["outdir"])
		else:
			plot_trans_arrays(conditions, transcript_arrays, args["outdir"])

		#Create web report
		FNULL = open(os.devnull, 'w')
		command = "unzip -o {0}/bootstrap-3.3.0-dist.zip -d {1}".format(path_to_stuff, args["outdir"]) #Move css and other stuff to output directory
		subprocess.call(command.split(), stdout=FNULL) 
		html = web_templates.create_transcript_html(transcripts) #Get HTML text 
		output = open(args["outdir"]+"/pynoncode.html", "w") #Write it out
		output.write(html)
		output.close()

		write_trans_summary(transcripts, args["outdir"])

	elif args["subparser_name"] == "frags":
		fragments = read_frag_input(args["input"], args["p"], args["pval"])
		transcripts = find_frag_transcripts(conditions, fragments, args["p"]) #Now contains coordinates of fragment, LFC and pvalue
		transcript_coords = preprocess_gtf(gtf, transcripts) #Annotation of transcripts
		transcript_arrays = read_directories_for_transcripts(conditions, transcript_coords, args["p"]) #Dict of numpy array containing counts of transcripts

		if args["a"]: #Average over conditions by reversing numpy array dict
			averaged_array = average_arrays(conditions, transcript_arrays, transcript_coords)
			plot_frag_arrays(conditions, averaged_array, args["outdir"], transcript_coords, transcripts, args["p"])
		else:
			plot_frag_arrays(conditions, transcript_arrays, args["outdir"], transcript_coords, transcripts, args["p"]) 

		#Creating web report
		FNULL = open(os.devnull, 'w')
		command = "unzip -o {0}/bootstrap-3.3.0-dist.zip -d {1}".format(path_to_stuff, args["outdir"]) #Move css and other stuff to output directory
		subprocess.call(command.split(), stdout=FNULL) 
		html = web_templates.create_fragment_html(transcripts, args["p"]) #Get HTML text 
		output = open(args["outdir"]+"/pynoncode.html", "w") #Write it out
		output.write(html)
		output.close()

		write_frag_summary(transcripts, args["p"], args["outdir"])
