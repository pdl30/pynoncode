#!/usr/bin/python

import subprocess
import re
import argparse


def stripper(fasta):
	result = {}
	with open(fasta) as f:
		for name, seq in fastq_fasta_parser.read_fasta(f):
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

def paired_bowtie(index, clipped=False):
	if clipped==False:
		sam1_o = open("unclipped_unique.sam", "wb")
		report1_o = open("unclipped_unique_report.txt", "wb")
		sam2_o = open("unclipped_multimap.sam", "wb")
		report2_o = open("unclipped_multi_report.txt", "wb")
		uniq = "bowtie --best -f -m 1 -v 2 --sam --un unclipped_unique_unmapped.fa {0} -1 original_fasta_1.fa -2 original_fasta_2.fa".format(index)
		multi= "bowtie --best -k 10 -f -m 500 -v 2 --sam --un unclipped_multi_unmapped.fa {0} -1 unclipped_unique_unmapped_1.fa -2 unclipped_unique_unmapped_2.fa".format(index)
		p = subprocess.Popen(uniq.split(), stdout = sam1_o, stderr=report1_o)
		p.communicate()
		p = subprocess.Popen(multi.split(), stdout = sam2_o, stderr=report2_o)
		p.communicate()
	else:
		sam1_o = open("clipped_unique.sam", "wb")
		report1_o = open("clipped_unique_report.txt", "wb")
		sam2_o = open("clipped_multimap.sam", "wb")
		report2_o = open("clipped_multimap_report.txt", "wb")
		uniq = "bowtie --best -f -m 1 -v 2 --sam --un clipped_unmapped.fa {0} -1 clipped_1.fa -2 clipped_2.fa".format(index)
		multi= "bowtie --best -k 10  -f -m 500 -v 2 --sam {0} -1 clipped_unmapped_1.fa -2 clipped_unmapped_2.fa".format(index)
		p = subprocess.Popen(uniq.split(), stdout = sam1_o, stderr=report1_o)
		p.communicate()
		p = subprocess.Popen(multi.split(), stdout = sam2_o, stderr=report2_o)
		p.communicate()

def single_bowtie(index, clipped=False):
	if clipped==False:
		sam1_o = open("unclipped_unique.sam", "wb")
		report1_o = open("unclipped_unique_report.txt", "wb")
		sam2_o = open("unclipped_multimap.sam", "wb")
		report2_o = open("unclipped_multi_report.txt", "wb")
		uniq = "bowtie --best -f -m 1 -v 2 --sam --un unclipped_unique_unmapped.fa {0} original_fasta.fa".format(index)
		multi= "bowtie --best -k 10 -f -m 500 -v 2 --sam --un unclipped_multi_unmapped.fa {0} unclipped_unique_unmapped.fa".format(index)
		p = subprocess.Popen(uniq.split(), stdout = sam1_o, stderr=report1_o)
		p.communicate()
		p = subprocess.Popen(multi.split(), stdout = sam2_o, stderr=report2_o)
		p.communicate()
	else:
		sam1_o = open("clipped_unique.sam", "wb")
		report1_o = open("clipped_unique_report.txt", "wb")
		sam2_o = open("clipped_multimap.sam", "wb")
		report2_o = open("clipped_multimap_report.txt", "wb")
		uniq = "bowtie --best -f -m 1 -v 2 --sam --un clipped_unique_unmapped.fa {0} clipped_fasta.fa".format(index)
		multi= "bowtie --best -k 10 -f -m 500 -v 2 --sam {0} clipped_unique_unmapped.fa".format(index)
		p = subprocess.Popen(uniq.split(), stdout = sam1_o, stderr=report1_o)
		p.communicate()
		p = subprocess.Popen(multi.split(), stdout = sam2_o, stderr=report2_o)
		p.communicate()

def grep_unique(samfile):
	out = re.sub(".sam", ".unique.sam", samfile)
	out2 = re.sub(".sam", ".multi.sam", samfile)
	output=  open(out, "w")
	output2=  open(out2, "w")
	with open(samfile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if line.startswith("@"):
				output.write("{}\n".format(line)),
				output2.write("{}\n".format(line)),
				continue
			if len(word) > 12:
				m = re.match("XS:i:", word[12])
				if m:
					if int(word[1]) == 147 or int(word[1]) == 83 or int(word[1]) == 99 or int(word[1]) == 163 or int(word[1]) == 81 or int(word[1]) == 97 or int(word[1]) == 145 or int(word[1]) == 161:
						output2.write("{}\n".format(line)),
				else:
					if int(word[1]) == 147 or int(word[1]) == 83 or int(word[1]) == 99 or int(word[1]) == 163 or int(word[1]) == 81 or int(word[1]) == 97 or int(word[1]) == 145 or int(word[1]) == 161:
						output.write("{}\n".format(line)),

def grep_single_unique(samfile):
	out = re.sub(".sam", ".unique.sam", samfile)
	out2 = re.sub(".sam", ".multi.sam", samfile)
	output=  open(out, "w")
	output2=  open(out2, "w")
	with open(samfile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if line.startswith("@"):
				output.write("{}\n".format(line)),
				output2.write("{}\n".format(line)),
				continue
			if len(word) > 12:
				m = re.match("XS:i:", word[12])
				if m:
					if int(word[1]) == 0 or int(word[1]) == 16:
						output2.write("{}\n".format(line)),
				else:
					if int(word[1]) == 0 or int(word[1]) == 16:
						output.write("{}\n".format(line)),


def paired_bowtie2(index, clipped=False):
	if clipped==False:
		report1_o = open("unclipped_unique_report.txt", "wb")
		uniq = "bowtie2 -k 10 -N 1 -f -p 12 --no-mixed --no-discordant --un-conc unmapped_round1.fa -x {0} -1 fasta_1.fa -2 fasta_2.fa -S tmp.sam".format(index)
		p = subprocess.Popen(uniq.split(), stderr=report1_o)
		p.communicate()
		grep_unique("tmp.sam")
		subprocess.call(["mv", "tmp.unique.sam", "bowtie2.uc.unique.sam"])
		subprocess.call(["mv", "tmp.multi.sam", "bowtie2.uc.multi.sam"])
	else:
		report1_o = open("clipped_unique_report.txt", "wb")
		uniq = "bowtie2 -k 10 -N 1 -f -p 12 --no-mixed --no-discordant --un-conc unmapped_round2.fa -x {0} -1 clipped_1.fa -2 clipped_2.fa -S tmp.sam".format(index)
		p = subprocess.Popen(uniq.split(), stderr=report1_o)
		p.communicate()
		grep_unique("tmp.sam")
		subprocess.call(["mv", "tmp.unique.sam", "bowtie2.c.unique.sam"])
		subprocess.call(["mv", "tmp.multi.sam", "bowtie2.c.multi.sam"])

def single_bowtie2(index, clipped=False):
	if clipped==False:
		report1_o = open("unclipped_unique_report.txt", "wb")
		uniq = "bowtie2 -k 10 -N 1 -f -p 12 --un unmapped_round1.fa -x {0} -U fasta.fa -S tmp.sam".format(index)
		p = subprocess.Popen(uniq.split(), stderr=report1_o)
		p.communicate()
		grep_single_unique("tmp.sam")
		subprocess.call(["mv", "tmp.unique.sam", "bowtie2.uc.unique.sam"])
		subprocess.call(["mv", "tmp.multi.sam", "bowtie2.uc.multi.sam"])
	else:
		report1_o = open("clipped_unique_report.txt", "wb")
		uniq = "bowtie2 -k 10 -N 1 -f -p 12 --un unmapped_round2.fa -x {0} -U clipped_fasta.fa -S tmp.sam".format(index)
		p = subprocess.Popen(uniq.split(), stderr=report1_o)
		p.communicate()
		grep_single_unique("tmp.sam")
		subprocess.call(["mv", "tmp.unique.sam", "bowtie2.c.unique.sam"])
		subprocess.call(["mv", "tmp.multi.sam", "bowtie2.c.multi.sam"])


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Runs bowtie\n')
	parser.add_argument('-p','--paired', help='Options True/False, is sample paired end?', required=False)
	parser.add_argument('-i','--index', help='Bowtie index', required=True)
	parser.add_argument('-c','--clipped', help='Options True/False, has sample been clipped?', required=False)
	args = vars(parser.parse_args())
	index = args["index"]
	if args["paired"] == True:
		if args["clipped"] == True:
			paired_bowtie(index, True)
		else:
			paired_bowtie(index, False)
	else:
		if args["clipped"] == True:
			single_bowtie(index, True)
		else:
			single_bowtie(index, False)
