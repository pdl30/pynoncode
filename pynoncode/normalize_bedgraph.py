import re

def normalize_bedgraph(bedgraph_file):
	total = 0
	data = {}
	out = re.sub(".bedGraph", "_tmp.bedGraph", bedgraph_file)
	output = open(out, "w")
	with open(bedgraph_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			total += float(word[3])
			id1 = "{}\t{}\t{}".format(word[0], word[1], word[2])
			data[id1] = word[3]
	for key in sorted(data.keys()):
		k = key.split("\t")
		new_count = float(data[key])/total
		rpm  = new_count * 1000000
		output.write("{}\t{}\t{}\t{}\n".format(k[0], k[1], k[2], rpm)),
