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
import ConfigParser
import sqlite3 as lite
import HTSeq
import pkg_resources

def read_directory(idir):
	transcript_counts = {}
	fragment_counts = {}
	with open(idir+"/transcripts.txt") as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			transcript_counts[word[0]] = round(float(word[1]))
	with open(idir+"/fragments.txt") as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			fragment_counts[word[0]] = round(float(word[1])), word[2], word[3], word[4], word[5], word[6]
	return transcript_counts, fragment_counts

def create_db(conditions, gtf, paired, outdir):
	seen = {}	
	anno_data = annotate_transcript(gtf)
	try:
		seen1 = {}
		seen2 = {}
		count_id2 = 0
		count_id1 = 0
		con = lite.connect(outdir +'/ncdiff.db')
		cur = con.cursor()  
		cur.execute("DROP TABLE IF EXISTS Transcripts;")
		cur.execute("DROP TABLE IF EXISTS Fragments;")
			#Create your new table
		cur.execute("CREATE TABLE Transcripts(Transcript_ID TEXT PRIMARY KEY, Chromosome TEXT, Start INT, End INT, Strand TEXT, biotype TEXT, Pvalue REAL, LFC REAL);")
		cur.execute("CREATE TABLE Fragments(Id INTEGER PRIMARY KEY, Fragment TEXT, Gene_ID TEXT, Chromosome TEXT, Start INT, End INT, Strand TEXT, biotype TEXT, Pvalue REAL, LFC REAL);")
		for sample in sorted(conditions.keys()):
			transcript_counts, fragment_counts = read_directory(sample)
			print("\nAdding {} fragments to database\n".format(sample)),
			#Add columns to your database
			cur.execute("alter table Transcripts add column '%s' 'int'" % sample)	
			cur.execute("alter table Fragments add column '%s' 'int'" % sample)
			for transcript in transcript_counts:
				#If fragment already in database, update the counts for that sample
				if transcript in seen2:
					tup = (transcript_counts[transcript], transcript, )
					cur.execute("UPDATE Transcripts SET '%s'=? WHERE Transcript_ID=?" % sample, (tup))
				#Otherwise, add a new row with the fragment information
				else:
					tup = (transcript, int(transcript_counts[transcript]), )
					cur.execute("INSERT INTO Transcripts(Transcript_ID, '%s') VALUES (?, ?)" % sample, (tup))
				seen2[transcript] = count_id2
				count_id2 += 1	

			for frag in fragment_counts:
				#If fragment already in database, update the counts for that sample
				if frag in seen1:
					tup = (fragment_counts[frag][1], seen1[frag], )
					cur.execute("UPDATE Fragments SET '%s'=? WHERE Id=?" % sample, (tup))
				#Otherwise, add a new row with the fragment information
				else:
				#	print fragment_counts[frag]
					cur.execute("INSERT INTO Fragments(Id, Fragment, Gene_ID, Chromosome, Start, End, Strand, '%s') VALUES (?, ?, ?, ?, ?, ?, ?, ?)" % sample, 
						(count_id1, frag, fragment_counts[frag][5], fragment_counts[frag][1],fragment_counts[frag][2], fragment_counts[frag][3], fragment_counts[frag][4], int(fragment_counts[frag][0])))
					seen1[frag] = count_id1
					count_id1 += 1	
			con.commit()
		for key in anno_data:
			tup = (anno_data[key][0], anno_data[key][1], anno_data[key][2], anno_data[key][3], anno_data[key][4], key, )
			cur.execute("UPDATE Transcripts SET Chromosome=?, Start=?, End=?, Strand=?, biotype=? WHERE Transcript_ID=?", tup)
	except lite.Error, e:
		if con:
			con.rollback()
			print "Error %s:" % e.args[0]
			sys.exit(1)
	finally:
		if con:
			con.close() 

def annotate_transcript(gtf):
	data = {}
	gtf_file = HTSeq.GFF_Reader(gtf)
	for feature in gtf_file:
		transcript = feature.attr.get('transcript_id', None)
		if transcript:
			pos, strand = str(feature.iv).split("/")
			chr1, pos1 = pos.split(":")
			pos1 = pos1[1:]
			pos1 = pos1[:-1]
			start, end = pos1.split(",")
			data[transcript] = chr1.rstrip(), start.rstrip(), end.rstrip(), strand.rstrip(), feature.source.rstrip()
	return data

def r_for_diff(table, paired, sample_dict, cond1, cond2, outdir, cutoff=0.1): ##NEEDS TO BE CHANGED! FOR TRANSCRIPT OUTPUT
	#Do I convert this to RPY2 or keep it as is?
	rscript =  "library(DESeq)\n"
	rscript += "library(RSQLite)\n"
	rscript += "library(biomaRt)\n"
	#Get data from database 
	rscript += "drv <- dbDriver('SQLite')\n"
	rscript += "con <- dbConnect(drv, '{}/ncdiff.db')\n".format(outdir) # Include this as variable? Not necessary at current time
	if table == "Transcripts":
		rscript += "full_table = dbGetQuery(con, 'Select * FROM Transcripts')\n"
	elif table == "Fragments":
		rscript += "full_table = dbGetQuery(con, 'Select * FROM Fragments')\n"

	rscript += "rownames(full_table) <- full_table[,1]\n"

	rscript += "p <- read.table('tmp_design.txt', header=F)\n"
	count = 0
	for key in sorted(sample_dict):
		if count == 0:
			rscript += "counts_data <- as.matrix(full_table${0}); colnames(counts_data) <- \"{0}\"\n".format(key)
		else:
			rscript += "counts_data <- cbind(counts_data, {0}=full_table${0})\n".format(key)
		count+=1
	rscript += "rownames(counts_data) <- full_table[,1]\n"
	rscript += "counts_data[is.na(counts_data)] <- 0\n" #Convert Nones to 0's
	rscript += "a <- rowSums(counts_data)\n"

	if table == "Transcripts":
		rscript += "filtered_counts <- counts_data[a>20,]\n" #??? what filter should I use?
	elif table == "Fragments":
		rscript += "filtered_counts <- counts_data[a>10,]\n"
	rscript += "pdata <- factor(p[,2])\n"
	#Initalise DESEQ with counts and pdata
	rscript += "cds <- newCountDataSet(filtered_counts,pdata)\n"
	rscript += "cds <- estimateSizeFactors(cds)\n"
	rscript += "set <- counts(cds,normalized=TRUE)\n"
	rscript += "cds <- estimateDispersions(cds,  fitType ='local')\n"
	rscript += "res <- nbinomTest(cds, '{0}', '{1}' )\n".format(cond1, cond2) #Run differental test
#	rscript += "resSig <- res[ res$padj < 0.1, ]\n" #Must include variable cutoff!
	rscript += "full_results <- data.frame(res, set)\n"
	rscript += "full_results$id <- NULL; full_results$baseMeanA <- NULL; full_results$baseMeanB <- NULL;\n" #To get rid of useless columns
	rscript += "id <- match(rownames(full_results), full_table[,1])\n"

	#Write to file
	if table == "Transcripts":
		rscript += "final_table <- cbind(ID=rownames(full_results), full_results)\n"
		rscript += "write.table(final_table, file='{0}vs{1}_transcript_analysis.csv', sep='\\t', quote=F, row.names=F)\n".format(cond1, cond2)
	elif table == "Fragments":
		if paired==True:
			rscript += "final_table <- cbind(ID=rownames(full_results), full_results, Pair1=full_table[id,2], Pair2=full_table[id,3], Mapped_gene=full_table[id,4])\n"
		else:
			rscript += "final_table <- cbind(ID=rownames(full_results), full_results, Fragment=full_table[id,2], Mapped_gene=full_table[id,3])\n"
		rscript += "write.table(final_table, file='{0}vs{1}_fragment_analysis.csv', sep='\\t', quote=F, row.names=F)\n".format(cond1, cond2)
	return rscript

def add_results_to_db(deseq_output, table):
	try:
		con = lite.connect('Fragments.db')
		cur = con.cursor()  
		print("\nAdding differental expression results to database!\n")
		with open(deseq_output) as f:
			next(f)
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				tup = (word[4], word[3], word[0], )
				if table == "Transcripts":
					cur.execute("UPDATE Transcripts SET Pvalue=?, LFC=? WHERE Transcript_ID=?", (tup))
				elif table== "Fragments":
					cur.execute("UPDATE Fragments SET Pvalue=?, LFC=? WHERE Id=?", (tup))
			con.commit()
	except lite.Error, e:
		if con:
			con.rollback()
			print "Error %s:" % e.args[0]
			sys.exit(1)
	finally:
		if con:
			con.close() 

def create_design_for_R(idict):
	output = open("tmp_design.txt", "w")
	for key in sorted(idict.keys()):
		output.write("{}\t{}\n".format(key, idict[key]))

def ConfigSectionMap(section):
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

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Processes ncRNA counts.\n')
	parser.add_argument('-config', help='Config file containing parameters, please see documentation for usage!', required=True)
	parser.add_argument('-gtf', help='GTF file, if not supplied, will use default', required=False)
	parser.add_argument('-paired', action='store_true', help='Are samples paired end?', required=False)
	parser.add_argument('-df', action='store_true', help='Perform differential expression analysis on uniquely aligned fragments', required=False)
	parser.add_argument('-out', help='Output results to this directory', required=False) #Will implement this later!
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	if args["gtf"]:
		gtf = args["gtf"]
	else:
		gtf = pkg_resources.resource_filename('ncpipe', 'data/mm10_ncRNA.gtf')
	
	conditions = ConfigSectionMap("Conditions")
	comparisons = ConfigSectionMap("Comparisons") #Can only support one comparison at a time! Doesn't work otherwise with the database!
	
	if not os.path.isdir(args["out"]):
		os.mkdir(args["out"])

	for comp in comparisons:
		c = comparisons[comp].split(",")
		comps = [x.strip(' ') for x in c]
	#	create_db(conditions, gtf, args["paired"], args["out"])
		create_design_for_R(conditions)
		rscript = r_for_diff("Transcripts", args["paired"], conditions, comps[0], comps[1], args["out"], cutoff=0.1)
		rcode = open(args["out"]+"/transcript_rcode.R", "w")
		rcode.write(rscript)
		rcode.close()