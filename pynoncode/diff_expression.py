#!/usr/bin/python

########################################################################
# 19 May 2014
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

def join_trans_counts(idict, outdir):
	print "==> Combining transcript counts...\n"
	transcripts = {}
	output = open(outdir+"/combined_transcript_counts.tsv", "w")
	output.write("ID"), #Sort out header
	for idir in sorted(idict):
		output.write("\t{}".format(idir)),
		with open(idir + "/transcript_counts.txt") as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				if word[0].startswith("__"):
					pass
				else:
					if word[0] in transcripts:
						transcripts[word[0]][idir] = word[1]
					else:
						transcripts[word[0]] = {}
						transcripts[word[0]][idir] = word[1]
	output.write("\n"),
	#First loop over transcripts
	for trans in sorted(transcripts):
		#Then loop over samples:
		output.write("{}".format(trans)),
		for idir in sorted(idict):
			value = round(float(transcripts[trans].get(idir, 0)))
			output.write("\t{}".format(value)),
		output.write("\n"),
	output.close()

def join_frag_counts(idict, outdir, paired):
	print "==> Combining fragment counts...\n"
	data = {}
	output = open(outdir+"/combined_fragment_counts.tsv", "w")
	output.write("ID"),
	for idir in sorted(idict):
		output.write("\t{}".format(idir)),
		if paired:
			with open(idir + "/fragment_counts.txt") as f:
				for line in f:
					line = line.rstrip()
					word = line.split("\t")
					if int(word[7]) == 1:
						next_line = next(f).rstrip()
						next_word = next_line.split("\t")
						if int(next_word[7]) == 2:
							if (word[3], next_word[3]) in data:
								if idir in data[(word[3], next_word[3])]:
									data[(word[3], next_word[3])][idir] += float(word[4]) #Allow for mulitple mapped fragments
								else:
									data[(word[3], next_word[3])][idir] = float(word[4])
							else:
								data[(word[3], next_word[3])][idir] = {} #New fragment
								data[(word[3], next_word[3])][idir] = float(word[4])
		else:
			with open(idir + "/fragment_counts.txt") as f:
				for line in f:
					line = line.rstrip()
					word = line.split("\t")
					if word[3] in data:
						if idir in data[word[3]]:
							data[word[3]][idir] += float(word[4])
						else:
							data[word[3]][idir] = float(word[4])
					else:
						data[word[3]] = {}
						data[word[3]][idir] = float(word[4])
	output.write("\n"),
	for frag in sorted(data):
		#Then loop over samples:
		if paired:
			output.write("{}\t{}".format(frag[0], frag[1])),
		else:
			output.write("{}".format(frag)),
		for idir in sorted(idict):
			value = round(float(data[frag].get(idir, 0)))
			output.write("\t{}".format(value)),
		output.write("\n"),
	output.close()


def write_deseq(sample_dict, cond1, cond2, padj, atype):
	print "==> Running differental expression analysis...\n"
	rscript =  "suppressMessages(library(DESeq2))\n"
	rscript += "pdata <- read.table('tmp_design.txt', header=T)\n"
	#Make sure they match!
	if atype == "trans":
		rscript += "counts <- read.table('combined_transcript_counts.tsv', sep='\\t', header=T, row.names=1)\n"
	elif atype == "frags":
		rscript += "counts <- read.table('combined_fragment_counts.tsv', sep='\\t', header=T, row.names=1)\n"
		rscript += "counts = counts[which(rowSums(counts) >= 10),]\n" #Hard cut-off, not sure??
	rscript += "rnaseq_dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(pdata), design = ~ condition)\n"
	rscript += "rnaseq_dds$condition <- factor(rnaseq_dds$condition, levels=unique(pdata[,2]))\n"
	rscript += "rnaseq_dds <- DESeq(rnaseq_dds)\n"
	rscript += "rnaseq_res <- results(rnaseq_dds, contrast=c('condition','{0}','{1}'))\n".format(cond1, cond2)
	rscript += "rnaseq_sig <- rnaseq_res[which(rnaseq_res$pvalue <= {}),]\n".format(padj)
	if atype == "trans":
		rscript += "write.table(rnaseq_sig, file='{0}_vs_{1}_sig_transcripts.tsv', sep='\\t', quote=F)\n".format(cond1, cond2)
		rscript += "write.table(rnaseq_res, file='{0}_vs_{1}_all_transcripts.tsv', sep='\\t', quote=F)\n".format(cond1, cond2)
	elif atype == "frags":
		rscript += "write.table(rnaseq_sig, file='{0}_vs_{1}_sig_fragments.tsv', sep='\\t', quote=F)\n".format(cond1, cond2)
		rscript += "write.table(rnaseq_res, file='{0}_vs_{1}_all_fragments.tsv', sep='\\t', quote=F)\n".format(cond1, cond2)
	return rscript

def create_design_for_R(outdir, idict):
	output = open(outdir + "/tmp_design.txt", "w")
	output.write("sampleName\tcondition\n"),
	for key in sorted(idict.keys()):
		output.write("{}\t{}\n".format(key, idict[key]))
	output.close()

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

def run_rcode(rscript, name):
	rcode = open(name, "w")
	rcode.write(rscript)
	rcode.close()
	try:
		subprocess.call(['Rscript', name])
	except:
		error("Error in running {}\n".format(name))
		error("Error: %s\n" % str(sys.exc_info()[1]))
		error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
		os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
		traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
		sys.exit(1)


def main():
	parser = argparse.ArgumentParser(description='Differential expression for RNA-seq experiments. Runs DESEQ2 by default\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")
	transcript_parser = subparsers.add_parser('trans', help="Runs differental expression for transcripts")
	fragment_parser = subparsers.add_parser('frags', help="Runs differental expression for fragments")

	transcript_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for examples of configuration and usage!', required=True)
	transcript_parser.add_argument('-a','--padj', help='Option for DESEQ2, default=0.05', default=0.05, required=False)
	transcript_parser.add_argument('-o','--outdir', help='Output directory', required=True)

	fragment_parser.add_argument('-c','--config', help='Config file containing parameters, please see documentation for examples of configuration and usage!', required=True)
	fragment_parser.add_argument('-a','--padj', help='Option for DESEQ2, default=0.05', default=0.05, required=False)
	fragment_parser.add_argument('-p','--paired', help='Are samples paired end?', action="store_true", required=False)
	fragment_parser.add_argument('-o','--outdir', help='Output directory', required=True)
	args = vars(parser.parse_args())
	conditions = []
	sample_dict = {}

	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	if os.path.isdir(args["outdir"]):
		print "Output directory already exists, existing results may be overwritten!"
	else:
		os.mkdir(args["outdir"])

		#Read design matrix and create list of conditions and directories
	conditions = ConfigSectionMap("Conditions", Config)
	comparisons = ConfigSectionMap("Comparisons", Config)
	create_design_for_R(args["outdir"], conditions) #Create design matrix

	if args["subparser_name"] == "trans":
		join_trans_counts(conditions, args["outdir"])
		os.chdir(args["outdir"])
		for comp in comparisons:
			c = comparisons[comp].split(",") #Names must be exact match for this to work!
			comps = [x.strip(' ') for x in c]
			rscript = write_deseq(conditions, comps[0], comps[1], args["padj"], args["subparser_name"]) ##Needs changing!!!
			run_rcode(rscript, "deseq2_rcode.R")

	elif args["subparser_name"] == "frags":
		join_frag_counts(conditions, args["outdir"], args["paired"])
		os.chdir(args["outdir"])
		for comp in comparisons:
			c = comparisons[comp].split(",") #Names must be exact match for this to work!
			comps = [x.strip(' ') for x in c]
			rscript = write_deseq(conditions, comps[0], comps[1], args["padj"], args["subparser_name"]) ##Needs changing!!!
			run_rcode(rscript, "deseq2_rcode.R")
