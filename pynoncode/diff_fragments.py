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
import sys
import re
import os
from collections import defaultdict
import ConfigParser
try:
	import sqlite3 as lite
except ImportError:
	import pip
	pip.main(['install', '--user', 'sqlite3'])
	import sqlite3 as lite
import HTSeq

def revcomp(dna, reverse):
	bases = 'ATGCNTACGN'
	complement_dict = {bases[i]:bases[i+5] for i in range(5)}
	if reverse:
		dna = reversed(dna)
	result = [complement_dict[base] for base in dna]
	return ''.join(result)

#Don't care about transcript positions, will annotate later!
def read_uniquemapped_comp(ifile, counts_dict, paired=False):
	counts=  {}
	frag_positions = {}
	transcript_counts = {}
	mapped_gene_fragment = {}
	input_file = open(ifile, "r")
	f = input_file.readlines()
	for index, line in enumerate(f):
		if index <= len(f) - 2:
			line = line.rstrip()
			word = line.split("\t")
			if paired==False:
				id1 = re.sub("ID:", "", word[3])
				if word[5] == "+":
					frag = word[4]
				else:
					frag = revcomp(word[4], True)

				frag_positions[frag] = [word[0], word[1], word[2], word[5]]
				mapped_gene_fragment[frag] = word[6]
				if frag not in counts:
					counts[frag] = int(counts_dict[id1])
				else:
					counts[frag] += int(counts_dict[id1])
				m = re.match("ambiguous", word[6])
				if word[6] == "no_feature" or m:
					pass
				else:
					if word[6] not in transcript_counts:
						transcript_counts[word[6]] = int(counts_dict[id1])	
					else:
						transcript_counts[word[6]] += int(counts_dict[id1])
			elif paired==True:
				next_line = f[index+1].rstrip()
				next_word = next_line.split("\t")
				if next_word[3] == word[3]:
					id1 = re.sub("ID:", "", word[3])

					if int(word[6]) == 1:
						frags = (word[4], next_word[4])
					elif int(next_word[6]) == 1:
						frags = (next_word[4], word[4])
					frag_positions[frags] = [word[0], word[1], word[2], word[5]]
					mapped_gene_fragment[frags] = word[7]
					
					if frags not in counts:
						counts[frags] = int(counts_dict[id1])
					else:
						counts[frags] += int(counts_dict[id1])
					m = re.match("ambiguous", word[7])
					if word[7] == "no_feature" or m:
						pass
					else:
						if word[7] not in transcript_counts:
							transcript_counts[word[7]] = int(counts_dict[id1])	
						else:
							transcript_counts[word[7]] += int(counts_dict[id1])
	return counts, transcript_counts, frag_positions, mapped_gene_fragment

def read_multimapped_comp(ifile, counts_dict, paired=False):
	transcript_counts = defaultdict(list)
	frag_positions = defaultdict(list)
	counts ={}
	mapped_gene_fragment =defaultdict(list)
	input_file = open(ifile, "r")
	f = input_file.readlines()
	for index, line in enumerate(f):
		if index <= len(f) - 2:
			if paired==False:
				line = line.rstrip()
				word = line.split("\t")

				id1 = re.sub("ID:", "", word[3])
				if word[5] == "+":
					frag = word[4]
				else:
					frag = revcomp(word[4], True)

				frag_positions[id1, word[6]].append([word[0], word[1], word[2], word[5], word[6]])
				mapped_gene_fragment[id1].append((frag, word[6]))

				m = re.match("ambiguous", word[6])
				if word[6] == "no_feature" or m:
					pass
				else:
					transcript_counts[id1].append((word[6], int(counts_dict[id1])))
			elif paired==True:
				line = line.rstrip()
				word = line.split("\t")
				next_line = f[index+1].rstrip()
				next_word = next_line.split("\t")
				id1 = re.sub("ID:", "", word[3])
				if int(word[6]) == 1:
					frags = (word[4], next_word[4])
				elif int(next_word[6]) == 1:
					frags = (next_word[4], word[4])

				frag_positions[id1, word[7]].append([word[0], word[1], word[2], word[5], word[7]])
				mapped_gene_fragment[id1].append((frags, word[7]))

				if next_word[3] == word[3] and word[7] == next_word[7]:
					m = re.match("ambiguous", word[6])
					if word[6] == "no_feature" or m:
						pass
					else:
						transcript_counts[id1].append((word[7], int(counts_dict[id1])))
	return transcript_counts, frag_positions, mapped_gene_fragment

def fragment_counts(idict, paired):
	#loops through the count files in listed directories and returns a dictionary of these counts and their fragments mapped genes
	#Could parallelise this but separate it from the create_db function. Would have to use the manager() in MP!

	unique_fragment_counts = {} #Counts of fragment
	transcript_counts = {} #Counts of transcripts
	fragment_positions = {}
	fragment_gene = {}
	multi_fragment_info = {}
	for name, sample1 in idict.items():
		used_multi_reads = open(sample1+"/used_multimapped_reads.txt", "w")
		read_counts = {}
		sample1 = sample1.rstrip()
		sample1 = re.sub("\"", "", sample1)
		unique = sample1+ "/unique_mapped.comp"
		multi = sample1+ "/multi_mapped.comp"
		counts_file = sample1+"/count_dict.txt"
		with open(counts_file) as f:
			for line in f:
				line = line.rstrip()
				word = line.split("\t")
				read_counts[word[0]] = word[1]

		unique_fragment_counts[name] = {} #Initialise sample dict
		transcript_counts[name] = {}
		fragment_positions[name] = {}
		fragment_gene[name] = {}
		multi_fragment_info[name] = {}
		multiple_fragment_newcounts = {}
		frag_total = {}
		if paired:
			unique_tmp, uniq_tmp_transcript_counts, uniq_pos_tmp, mapped_gene_fragment = read_uniquemapped_comp(unique, read_counts, paired=True)
			unique_fragment_counts[name].update(unique_tmp)
			multi_tmp, multi_fragment_positions, multi_fragment_gene = read_multimapped_comp(multi, read_counts, paired=True)
		else:
			unique_tmp, uniq_tmp_transcript_counts, uniq_pos_tmp, mapped_gene_fragment = read_uniquemapped_comp(unique, read_counts, paired=False)
			unique_fragment_counts[name].update(unique_tmp)
			multi_tmp, multi_fragment_positions, multi_fragment_gene = read_multimapped_comp(multi, read_counts, paired=False)

		read_total = {}
		multi_count = {}
		for read_id in multi_tmp:
			for trans, read_count in multi_tmp[read_id]:
				if trans not in uniq_tmp_transcript_counts:
					pass
				else:
					if read_id not in read_total:
						read_total[read_id] = int(uniq_tmp_transcript_counts[trans])
						multi_count[read_id] = 1
					else:
						read_total[read_id] += int(uniq_tmp_transcript_counts[trans])
						multi_count[read_id] += 1

		for read_id in multi_fragment_gene:
			for frag, trans in multi_fragment_gene[read_id]:
				if trans not in uniq_tmp_transcript_counts:
					pass
				else:
					if read_id not in frag_total:
						frag_total[read_id] = int(uniq_tmp_transcript_counts[trans])
					else:
						frag_total[read_id] += int(uniq_tmp_transcript_counts[trans])
		
		uniq_tmp_transcript_counts_copy = uniq_tmp_transcript_counts.copy()

		for read_id in multi_tmp:
			for trans, read_count in multi_tmp[read_id]:
				if trans not in uniq_tmp_transcript_counts:
					pass
				elif read_total[read_id] < 10:
					tmp_count3 = float(read_count)/float(multi_count[read_id])
					uniq_tmp_transcript_counts_copy[trans] += tmp_count3
					used_multi_reads.write("{}\t{}\t{}\n".format(read_id, trans, tmp_count3)),
				else:
					tmp_count = float(uniq_tmp_transcript_counts[trans])/read_total[read_id]
					tmp_count2 = tmp_count*read_count
					uniq_tmp_transcript_counts_copy[trans] += int(tmp_count2)
					used_multi_reads.write("{}\t{}\t{}\n".format(read_id, trans, tmp_count2)),

		tmp_c = 0
		for read_id in multi_fragment_gene:
			for frag, trans in multi_fragment_gene[read_id]:
				if trans not in uniq_tmp_transcript_counts:
					pass
				else:
					#print name, read_id, trans, uniq_tmp_transcript_counts[trans], frag_total[read_id], read_counts[read_id]
					tmp_count = float(uniq_tmp_transcript_counts[trans])/int(frag_total[read_id])
					tmp_count2 = tmp_count*int(read_counts[read_id])
					multiple_fragment_newcounts[tmp_c] = ((frag, tmp_count2, multi_fragment_positions[read_id, trans]))
					tmp_c += 1

		fragment_gene[name].update(mapped_gene_fragment)
		transcript_counts[name].update(uniq_tmp_transcript_counts_copy)
		fragment_positions[name].update(uniq_pos_tmp)
		multi_fragment_info[name].update(multiple_fragment_newcounts)

		del read_counts, uniq_tmp_transcript_counts, multi_tmp, read_total, uniq_pos_tmp, mapped_gene_fragment, uniq_tmp_transcript_counts_copy, multiple_fragment_newcounts
		used_multi_reads.close()
	return unique_fragment_counts, transcript_counts, fragment_positions, fragment_gene, multi_fragment_info

def create_db(idict, paired):
	if paired == True:
		unique_fragment_counts, transcript_counts, fragment_positions, fragment_gene, multi_fragment_info = fragment_counts(idict, True)
	else:
		unique_fragment_counts, transcript_counts, fragment_positions, fragment_gene, multi_fragment_info = fragment_counts(idict, False)
	seen = {}	
	try:
		seen1 = {}
		seen2 = {}
		seen3 = {}
		count_id2 = 0
		count_id1 = 0
		count_id3 = 0
		con = lite.connect('Fragments.db')
		cur = con.cursor()  
		cur.execute("DROP TABLE IF EXISTS Transcript_Counts;")

		cur.execute("DROP TABLE IF EXISTS Unique_Fragment_Counts;")
			#Create your new table
		cur.execute("CREATE TABLE Transcript_Counts(Transcript_ID TEXT PRIMARY KEY, Chromosome TEXT, Start INT, End INT, Strand TEXT, biotype TEXT, Pvalue REAL, LFC REAL);")
		if paired == True:
			cur.execute("CREATE TABLE Unique_Fragment_Counts(Id INTEGER PRIMARY KEY, Frag1 TEXT, Frag2 TEXT, Gene_ID TEXT, Chromosome TEXT, Start INT, End INT, Strand TEXT, biotype TEXT, Pvalue REAL, LFC REAL);")
		else:
			cur.execute("CREATE TABLE Unique_Fragment_Counts(Id INTEGER PRIMARY KEY, Frag1 TEXT, Gene_ID TEXT, Chromosome TEXT, Start INT, End INT, Strand TEXT, biotype TEXT, Pvalue REAL, LFC REAL);")
		for sample in sorted(idict.keys()):
	
			print("\nAdding {} fragments to database\n".format(sample)),
			#Add columns to your database
			cur.execute("alter table Transcript_Counts add column '%s' 'int'" % sample)	
			cur.execute("alter table Unique_Fragment_Counts add column '%s' 'int'" % sample)
			for transcript in transcript_counts[sample]:
				#If fragment already in database, update the counts for that sample
				if transcript in seen2:
					tup = (transcript_counts[sample][transcript], transcript, )
					cur.execute("UPDATE Transcript_Counts SET '%s'=? WHERE Transcript_ID=?" % sample, (tup))
				#Otherwise, add a new row with the fragment information
				else:
					tup = (transcript, int(transcript_counts[sample][transcript]), )
					cur.execute("INSERT INTO Transcript_Counts(Transcript_ID, '%s') VALUES (?, ?)" % sample, (tup))
				seen2[transcript] = count_id2
				count_id2 += 1	

			for frag in unique_fragment_counts[sample]:
				#If fragment already in database, update the counts for that sample
				if frag in seen1:
					tup = (unique_fragment_counts[sample][frag], seen1[frag], )
					cur.execute("UPDATE Unique_Fragment_Counts SET '%s'=? WHERE Id=?" % sample, (tup))
				#Otherwise, add a new row with the fragment information
				else:
					if paired == True:
						frag1, frag2 = frag
						coords = fragment_positions[sample][frag]
						cur.execute("INSERT INTO Unique_Fragment_Counts(Id, Frag1, Frag2, Gene_ID, Chromosome, Start, End, Strand, '%s') VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)" % sample, 
							(count_id1, frag1, frag2, fragment_gene[sample][frag], coords[0], coords[1], coords[2], coords[3], int(unique_fragment_counts[sample][frag])))
					else:
						coords = fragment_positions[sample][frag]
						cur.execute("INSERT INTO Unique_Fragment_Counts(Id, Frag1, Gene_ID, Chromosome, Start, End, Strand, '%s') VALUES (?, ?, ?, ?, ?, ?, ?, ?)" % sample, 
							(count_id1, frag, fragment_gene[sample][frag], coords[0], coords[1], coords[2], coords[3], int(unique_fragment_counts[sample][frag])))
					seen1[frag] = count_id1
					count_id1 += 1	

			for key in multi_fragment_info[sample]:
				frag, count, frag_pos = multi_fragment_info[sample][key]
				if frag in seen3:
					tup = ( int(count), seen3[frag], )
					cur.execute("UPDATE Unique_Fragment_Counts SET '%s'=? WHERE Id=?" % sample, (tup))
				else:
					if paired == True:
						frag1, frag2 = frag
						cur.execute("INSERT INTO Unique_Fragment_Counts(Id, Frag1, Frag2, Gene_ID, Chromosome, Start, End, Strand, '%s') VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)" % sample, 
							(count_id1, frag1, frag2,frag_pos[0][4], frag_pos[0][0], frag_pos[0][1], frag_pos[0][2], frag_pos[0][3], int(count)))
					else:
						cur.execute("INSERT INTO Unique_Fragment_Counts(Id, Frag1, Gene_ID, Chromosome, Start, End, Strand, '%s') VALUES (?, ?, ?, ?, ?, ?, ?, ?)" % sample, 
							(count_id1, frag, frag_pos[0][4], frag_pos[0][0], frag_pos[0][1], frag_pos[0][2], frag_pos[0][3], int(count)))
					count_id1 += 1	
					seen3[frag] = count_id1
			con.commit()
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
	try:
		con = lite.connect('Fragments.db')
		cur = con.cursor()  

		for key in data:
			tup = (data[key][0], data[key][1], data[key][2], data[key][3], data[key][4], key, )
			cur.execute("UPDATE Transcript_Counts SET Chromosome=?, Start=?, End=?, Strand=?, biotype=? WHERE Transcript_ID=?", tup)
		con.commit()
	except lite.Error, e:
		if con:
			con.rollback()
			print "Error %s:" % e.args[0]
			sys.exit(1)
	finally:
		if con:
			con.close() 

def R_for_diff(table, paired, sample_dict, cond1, cond2, cutoff=0.1): ##NEEDS TO BE CHANGED! FOR TRANSCRIPT OUTPUT
	#Do I convert this to RPY2 or keep it as is?
	rscript =  "library(DESeq)\n"
	rscript += "library(RSQLite)\n"
	rscript += "library(biomaRt)\n"
	#Get data from database 
	rscript += "drv <- dbDriver('SQLite')\n"
	rscript += "con <- dbConnect(drv, 'Fragments.db')\n" # Include this as variable? Not necessary at current time
	if table == "Transcript_Counts":
		rscript += "full_table = dbGetQuery(con, 'Select * FROM Transcript_Counts')\n"
	elif table == "Unique_Fragment_Counts":
		rscript += "full_table = dbGetQuery(con, 'Select * FROM Unique_Fragment_Counts')\n"

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

	if table == "Transcript_Counts":
		rscript += "filtered_counts <- counts_data[a>20,]\n" #??? what filter should I use?
	elif table == "Unique_Fragment_Counts":
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
	if table == "Transcript_Counts":
		rscript += "final_table <- cbind(ID=rownames(full_results), full_results)\n"
		rscript += "write.table(final_table, file='{0}vs{1}_transcript_analysis.csv', sep='\\t', quote=F, row.names=F)\n".format(cond1, cond2)
	elif table == "Unique_Fragment_Counts":
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
				if table == "Transcript_Counts":
					cur.execute("UPDATE Transcript_Counts SET Pvalue=?, LFC=? WHERE Transcript_ID=?", (tup))
				elif table== "Unique_Fragment_Counts":
					cur.execute("UPDATE Unique_Fragment_Counts SET Pvalue=?, LFC=? WHERE Id=?", (tup))
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
	parser.add_argument('-c','--CONFIG', help='Config file containing parameters, please see documentation for usage!', required=True)
	parser.add_argument('-g','--GTF', help='GTF file.', required=True)
	parser.add_argument('-p', action='store_true', help='Are samples paired end?', required=False)
	parser.add_argument('-f', action='store_true', help='Perform differential expression analysis on uniquely aligned fragments', required=False)
	parser.add_argument('-o','--OUTDIR', help='Output results to this directory, don\'t use! not implemented yet', required=False) #Will implement this later!
	args = vars(parser.parse_args())
	conditions = []
	sample_dict = {}
	count_id = 0
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["CONFIG"])

	#Read design matrix and create list of conditions and directories
	sample_dict = ConfigSectionMap("Sample_Names")
	conditions = ConfigSectionMap("Sample_Conditions")
	comparisons = ConfigSectionMap("Comparisons") #Can only support one comparison at a time! Doesn't work otherwise with the database!

	
	#Write R code for DESEQ analysis
	for comp in comparisons:
		c = comparisons[comp].split(",")
		comps = [x.strip(' ') for x in c]
		if args["p"]:
			create_db(sample_dict, True)
			create_design_for_R(conditions)
			annotate_transcript(args["GTF"])
			rscript = R_for_diff("Transcript_Counts", True, conditions, comps[0], comps[1]) ##Needs changing!!!
			rcode = open("transcript_rcode.R", "w")
			rcode.write(rscript)
			rcode.close()
			try:
				sts = subprocess.call(['Rscript', 'transcript_rcode.R'])
			except:
				error("Error in running differential analysis!\n")
				error("Error: %s\n" % str(sys.exc_info()[1]))
				error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
				os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
				traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
				sys.exit(1)

			#Add results to database
			
			results_file = "{0}vs{1}_transcript_analysis.csv".format(comps[0], comps[1])
			add_results_to_db(results_file, "Transcript_Counts")

			if args["f"]:
				rscript2 = R_for_diff("Unique_Fragment_Counts", True, conditions, comps[0], comps[1])
				rcode = open("fragment_rcode.R", "w")
				rcode.write(rscript2)
				rcode.close()
				try:
					sts = subprocess.call(['Rscript', 'fragment_rcode.R'])
				except:
					error("Error in running differential analysis!\n")
					error("Error: %s\n" % str(sys.exc_info()[1]))
					error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
					os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
					traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
					sys.exit(1)
				results_file = "{0}vs{1}_fragment_analysis.csv".format(comps[0], comps[1])
				add_results_to_db(results_file, "Unique_Fragment_Counts")

		else:
			create_db(sample_dict, False)
			create_design_for_R(conditions)
			annotate_transcript(args["GTF"])
			rscript = R_for_diff("Transcript_Counts", False, conditions, comps[0], comps[1])
			rcode = open("transcript_rcode.R", "w")
			rcode.write(rscript)
			rcode.close()
			try:
				sts = subprocess.call(['Rscript', 'transcript_rcode.R'])
			except:
				error("Error in running differential analysis!\n")
				error("Error: %s\n" % str(sys.exc_info()[1]))
				error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
				os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
				traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
				sys.exit(1)

			#Add results to database
			results_file = "{0}vs{1}_transcript_analysis.csv".format(comps[0], comps[1])
			add_results_to_db(results_file, "Transcript_Counts")

			if args["f"]:
				rscript2 = R_for_diff("Unique_Fragment_Counts", False, conditions, comps[0], comps[1])
				rcode = open("fragment_rcode.R", "w")
				rcode.write(rscript2)
				rcode.close()
				try:
					sts = subprocess.call(['Rscript', 'fragment_rcode.R'])
				except:
					error("Error in running differential analysis!\n")
					error("Error: %s\n" % str(sys.exc_info()[1]))
					error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
					os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
					traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
					sys.exit(1)
				results_file = "{0}vs{1}_fragment_analysis.csv".format(comps[0], comps[1])
				add_results_to_db(results_file, "Unique_Fragment_Counts")