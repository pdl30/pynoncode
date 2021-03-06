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

def revcomp(dna, reverse):
	bases = 'ATGCNTACGN'
	complement_dict = {bases[i]:bases[i+5] for i in range(5)}
	if reverse:
		dna = reversed(dna)
	result = [complement_dict[base] for base in dna]
	return ''.join(result)

def single_uniquemapped_comp(ifile, counts_dict):
	results = {}
	results["transcript_counts"] = {}
	results["fragment_counts"] = {}
	results["fragment_positions"] = {}
	results["transcript_id"] = {}
	input_file = open(ifile, "r")
	f = input_file.readlines()
	for index, line in enumerate(f):
		line = line.rstrip()
		word = line.split("\t")
		id1 = re.sub("ID:", "", word[3])
		m = re.match("ambiguous*", word[6])
		if word[6] == "no_feature" or m:
			pass
		else:
			if word[5] == "+":
				frag = word[4]
			else:
				frag = revcomp(word[4], True)
			results["fragment_positions"][frag] = [word[0], word[1], word[2], word[5]]
			results["transcript_id"][frag] = word[6]
			if frag not in results["fragment_counts"]:
				results["fragment_counts"][frag] = int(counts_dict[id1])
			else:
				results["fragment_counts"][frag] += int(counts_dict[id1])
		
			if word[6] not in results["transcript_counts"]:
				results["transcript_counts"][word[6]] = int(counts_dict[id1])	
			else:
				results["transcript_counts"][word[6]] += int(counts_dict[id1])
	return results


def paired_uniquemapped_comp(ifile, counts_dict):
	results = {}
	results["transcript_counts"] = {}
	results["fragment_counts"] = {}
	results["fragment_positions"] = {}
	results["transcript_id"] = {}
	input_file = open(ifile, "r")
	f = input_file.readlines()
	len_f = len(f)
	for i in range(len_f-1): 
		line = f[i].rstrip()
		word = line.split("\t")
		if i < len_f:
			next_line = f[i+1].rstrip()
		else:
			next_line = ""
		next_word = next_line.split("\t")
		if next_word[3] == word[3]: #This makes sure counts are only added once per pair!
			id1 = re.sub("ID:", "", word[3])
			m = re.match("ambiguous", word[7])
			if word[7] == "no_feature" or m:
				pass
			else:
				if int(word[6]) == 1:
					frags = (word[4], next_word[4])
				elif int(next_word[6]) == 1:
					frags = (next_word[4], word[4])
				results["fragment_positions"][frags] = [word[0], word[1], word[2],word[5], next_word[0],next_word[1],next_word[2], next_word[5]]
				results["transcript_id"][frags] = word[7]
				if frags not in results["fragment_counts"]:
					results["fragment_counts"][frags] = int(counts_dict[id1])
				else:
					results["fragment_counts"][frags] += int(counts_dict[id1])
				
				if word[7] not in results["transcript_counts"]:
					results["transcript_counts"][word[7]] = int(counts_dict[id1])	
				else:
					results["transcript_counts"][word[7]] += int(counts_dict[id1])
	return results

def single_multimapped_comp(ifile, counts_dict):
	results = {}
	results["transcript_counts"] = defaultdict(list)
	results["fragment_positions"] = defaultdict(list)
	results["transcript_id"] = defaultdict(list)
	input_file = open(ifile, "r")
	f = input_file.readlines()
	for index, line in enumerate(f):
		line = line.rstrip()
		word = line.split("\t")
		#Key is the read ID
		id1 = re.sub("ID:", "", word[3])
		if word[5] == "+":
			frag = word[4]
		else:
			frag = revcomp(word[4], True)
		results["fragment_positions"][id1, word[6]].append([word[0], word[1], word[2], word[5], word[6]])
		results["transcript_id"][id1].append((frag, word[6]))
		m = re.match("ambiguous", word[6])
		if word[6] == "no_feature" or m:
			pass
		else:
			results["transcript_counts"][id1].append((word[6], int(counts_dict[id1]))) #Tuple of the transcript id and read count
	return results

def paired_multimapped_comp(ifile, counts_dict):
	results = {}
	results["transcript_counts"] = defaultdict(list)
	results["fragment_positions"] = defaultdict(list)
	results["transcript_id"] = defaultdict(list)
	input_file = open(ifile, "r")
	f = input_file.readlines()
	len_f = len(f)
	for i in range(len_f-1):
		line = f[i].rstrip()
		word = line.split("\t")
		if i < len_f:
			next_line = f[i+1].rstrip()
		else:
			next_line = ""
		next_word = next_line.split("\t")
		id1 = re.sub("ID:", "", word[3])
		if next_word[3] == word[3] and word[7] == next_word[7]:
			m = re.match("ambiguous", word[6])
			if word[6] == "no_feature" or m:
				pass
			else:
				results["transcript_counts"][id1].append((word[7], int(counts_dict[id1])))
			if int(word[6]) == 1:
				frags = (word[4], next_word[4])
			elif int(next_word[6]) == 1:
				frags = (next_word[4], word[4])
			results["fragment_positions"][id1, word[7]].append([word[0], word[1], word[2], word[5], next_word[0],next_word[1],next_word[2], next_word[5]])
			results["transcript_id"][id1].append((frags, word[7]))
	return results

def read_count_file():
	counts_file = "count_dict.txt"
	read_counts = {}
	with open(counts_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			read_counts[word[0]] = word[1]
	return read_counts

def get_unique_total_for_reads(unique_dict, multi_dict):
	read_total = {}
	multi_count = {}
	#Gets total unique count for multimapped reads and multi count tells me how many genes are found in the unique dict
	for read_id in multi_dict["transcript_counts"]:
		for value in multi_dict["transcript_counts"][read_id]: #Transcript_counts is a tuple of transript and read count
			if value[0] not in unique_dict["transcript_counts"]: #only need transcript
				pass
			else:
				if read_id not in read_total:
					read_total[read_id] = int(unique_dict["transcript_counts"][value[0]])
					multi_count[read_id] = 1
				else:
					read_total[read_id] += int(unique_dict["transcript_counts"][value[0]])
					multi_count[read_id] += 1
	return read_total, multi_count

def distribute_transcripts(unique_dict, multi_dict, new_transcript_counts, unique_totals, multi_uniques_count):
	used_multi_reads = open("used_multimapped_reads.txt", "w")
	for read_id in multi_dict["transcript_counts"]:
		for trans, read_count in multi_dict["transcript_counts"][read_id]:
			if trans not in unique_dict["transcript_counts"]:
				pass
			elif unique_totals[read_id] < 10: #If the unique count for the gene is less than 10, distribute the multi-mapped equally between all the transcripts
				tmp_count3 = float(read_count)/float(multi_uniques_count[read_id])
				new_transcript_counts[trans] += tmp_count3
				used_multi_reads.write("{}\t{}\t{}\n".format(read_id, trans, tmp_count3)),
			else:
				tmp_count = float(unique_dict["transcript_counts"][trans])/unique_totals[read_id]
				tmp_count2 = tmp_count*read_count
				new_transcript_counts[trans] += int(tmp_count2)
				used_multi_reads.write("{}\t{}\t{}\n".format(read_id, trans, tmp_count2)),
	return new_transcript_counts

def distribute_fragments(unique_dict, multi_dict, unique_totals, count_dict):
	tmp_c = 0
	distributed_fragments = {}
	for read_id in multi_dict["transcript_id"]:
		for frag, trans in multi_dict["transcript_id"][read_id]:
			if trans not in unique_dict["transcript_counts"]:
				pass
			else:
				tmp_count = float(unique_dict["transcript_counts"][trans])/int(unique_totals[read_id])
				tmp_count2 = tmp_count*int(count_dict[read_id])
				if tmp_count2 >= 1:
					distributed_fragments[tmp_c] = ((frag, tmp_count2, multi_dict["fragment_positions"][read_id, trans], trans))
					tmp_c += 1
	return distributed_fragments

def cleanup(orig_dir, input_dir):
	os.chdir(orig_dir)
	os.chdir(input_dir)
	command = "mv pynoncode/transcript_counts.txt ."
	command2 = "mv pynoncode/fragment_counts.txt ."
	subprocess.call(command.split())
	subprocess.call(command2.split())

def write_transcripts(distributed_transcripts):
	output = open("transcript_counts.txt", "w")
	for transcript in distributed_transcripts:
		output.write("{}\t{}\n".format(transcript, distributed_transcripts[transcript])),
	output.close()

def write_multiple_fragments(distributed_fragments, paired):
	output1 = open("fragment_counts.txt", "w")
	for fragment in distributed_fragments:
			#First Read
		if paired:
			output1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t1\n".format(distributed_fragments[fragment][2][0][0],distributed_fragments[fragment][2][0][1],distributed_fragments[fragment][2][0][2],
				distributed_fragments[fragment][0][0],distributed_fragments[fragment][1], distributed_fragments[fragment][2][0][3],distributed_fragments[fragment][3])),
		#Second Read
			output1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t2\n".format(distributed_fragments[fragment][2][0][4],distributed_fragments[fragment][2][0][5],distributed_fragments[fragment][2][0][6],
				distributed_fragments[fragment][0][1],distributed_fragments[fragment][1], distributed_fragments[fragment][2][0][7],distributed_fragments[fragment][3])),
		else:
			output1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(distributed_fragments[fragment][2][0][0],distributed_fragments[fragment][2][0][1],distributed_fragments[fragment][2][0][2],
				distributed_fragments[fragment][0],distributed_fragments[fragment][1], distributed_fragments[fragment][2][0][3],distributed_fragments[fragment][3])),
	output1.close()

def write_single_fragments(unique_data, paired):
	output1 = open("fragment_counts.txt", "a")
	for fragment in unique_data["fragment_counts"]:
		if paired:
			output1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t1\n".format(unique_data["fragment_positions"][fragment][0], unique_data["fragment_positions"][fragment][1], 
				unique_data["fragment_positions"][fragment][2], fragment[0], unique_data["fragment_counts"][fragment], unique_data["fragment_positions"][fragment][3],
				unique_data["transcript_id"][fragment])),
			output1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t2\n".format(unique_data["fragment_positions"][fragment][4], unique_data["fragment_positions"][fragment][5], 
				unique_data["fragment_positions"][fragment][6], fragment[1], unique_data["fragment_counts"][fragment], unique_data["fragment_positions"][fragment][7],
				unique_data["transcript_id"][fragment])),
		else:
			output1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(unique_data["fragment_positions"][fragment][0], unique_data["fragment_positions"][fragment][1], 
				unique_data["fragment_positions"][fragment][2], fragment, unique_data["fragment_counts"][fragment], unique_data["fragment_positions"][fragment][3],
				unique_data["transcript_id"][fragment])),
	output1.close()

def main():
	parser = argparse.ArgumentParser(description='Processes ncRNA samples from fastq files to sam file.\n Please ensure FASTQ files are in current directory.\n ')
	parser.add_argument('-i', '--input', help='Input directory after ncalign has been run', required=True)
	parser.add_argument('-p', help='Use if experiment is paired end', action="store_true", required=False)
	parser.add_argument('-m', action='store_true', help='Use multiple mapped reads and add to final count for genes', required=False)
	args = vars(parser.parse_args())
	orig_path = os.getcwd()
	path = os.path.join(args["input"], "pynoncode")
	os.chdir(path)
	read_counts = read_count_file()

	#Data dicts include transcript_counts, transcript_id, fragment_positions, fragment_counts [unique only!]
	if args["p"]:
		unique_data = paired_uniquemapped_comp("unique_mapped.BED", read_counts)
	else:
		unique_data = single_uniquemapped_comp("unique_mapped.BED", read_counts)

	if args["m"]:
		if args["p"]:
			multi_data = paired_multimapped_comp("multi_mapped.BED", read_counts)
		else:
			multi_data = single_multimapped_comp("multi_mapped.BED", read_counts)
		unique_totals, multi_uniques_count = get_unique_total_for_reads(unique_data, multi_data)
		
		new_transcript_counts = unique_data["transcript_counts"].copy() #Copy original dict for new data

		distributed_transcripts = distribute_transcripts(unique_data, multi_data, new_transcript_counts, unique_totals, multi_uniques_count)
		distributed_fragments = distribute_fragments(unique_data, multi_data, unique_totals, read_counts)

		write_transcripts(distributed_transcripts)
		write_multiple_fragments(distributed_fragments, args["p"])
		
	else:
		write_transcripts(unique_data["transcript_counts"])

	write_single_fragments(unique_data, paired)
	cleanup(orig_path, args["input"])