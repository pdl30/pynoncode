#!/usr/bin/python

########################################################################
# 31 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

from pynoncode import basic_template

def trans_table(transcript_info):
	thead = """
<thead>
<tr>
<th>Transcript Name</th>
<th>P Value</th>
<th>Log Fold Change</th>
<th>Plot of sample coverage</th>
</tr>
</thead>"""
	count = 0 
	for trans in sorted(transcript_info):
		link = "plots/{}.png".format(trans)
		if count == 0:
			body2 = """
<tbody>
<tr>
<td>{0}</td>
<td>{1}</td>
<td>{2}</td>
<td><a href="{3}">Link</a></td>
</tr>
        	""".format(trans, transcript_info[trans][1], transcript_info[trans][0], link)
		else:
			body2 += """
<tr>
<td>{0}</td>
<td>{1}</td>
<td>{2}</td>
<td><a href="{3}">Link</a></td>
</tr>
            """.format(trans, transcript_info[trans][1], transcript_info[trans][0], link)
		count+=1 
	table = thead
	table += body2
	return table


def frags_table(transcipts, paired):
	if paired:
		thead = """
<thead>
<tr>
<th>Chromosome</th>
<th>Start</th>
<th>End</th>
<th>Read1</th>
<th>Chromosome</th>
<th>Start</th>
<th>End</th>
<th>Read2</th>
<th>P Value</th>
<th>Log Fold Change</th>
<th>Mapped Transcript</th>
<th>Plot of sample coverage</th>
</tr>
</thead>
<tbody>"""
		for trans in sorted(transcipts):
			a = 1
			for frag_pos in sorted(transcipts[trans]): #Now a list containing fragment, lfc, pvalue, padj
				print trans, frag_pos
				link = "plots/{}_{}.png".format(trans, a)
				a += 1
				thead  += """

<tr>
<td>{0}</td>
<td>{1}</td>
<td>{2}</td>
<td>{3}</td>
<td>{4}</td>
<td>{5}</td>
<td>{6}</td>
<td>{7}</td>
<td>{8}</td>
<td>{9}</td>
<td>{10}</td>
<td><a href="{11}">Link</a></td>
</tr>
 	  			""".format(frag_pos[0], frag_pos[1], frag_pos[2], transcipts[trans][frag_pos][0], frag_pos[0], frag_pos[3], frag_pos[4], transcipts[trans][frag_pos][1], 
 	  				transcipts[trans][frag_pos][2][1], transcipts[trans][frag_pos][2][0], trans, link)
	else: #Unpaired samples
		thead = """
<thead>
<tr>
<th>Chromosome</th>
<th>Start</th>
<th>End</th>
<th>Read1</th>
<th>P Value</th>
<th>Log Fold Change</th>
<th>Mapped Transcript</th>
<th>Plot of sample coverage</th>
</tr>
</thead>
<tbody>"""
		for trans in sorted(transcipts):
			a = 1
			for frag_pos in sorted(transcipts[trans]): #Now a list containing fragment, lfc, pvalue, padj
				link = "plots/{}_{}.png".format(trans, a)
				a += 1
				thead  += """

<tr>
<td>{0}</td>
<td>{1}</td>
<td>{2}</td>
<td>{3}</td>
<td>{4}</td>
<td>{5}</td>
<td>{6}</td>
<td><a href="{7}">Link</a></td>
</tr>
 	  			""".format(frag_pos[0], frag_pos[1], frag_pos[2], transcipts[trans][frag_pos][0], transcipts[trans][frag_pos][1][1], transcipts[trans][frag_pos][1][0], trans, link)
	return thead



def create_transcript_html(transcripts):
	total_string = basic_template.t_header()
	total_string += trans_table(transcripts)
	total_string += basic_template.ending()
	return total_string

def create_fragment_html(transcripts, paired):
	total_string = basic_template.f_header()
	total_string += frags_table(transcripts, paired)
	total_string += basic_template.ending()
	return total_string


def counts_trans_table(transcript_counts):
	c = 0
	samples = []
	for trans in transcript_counts:
		for sample in transcript_counts[trans]:
			if c == 0:
				samples.append(sample)
				c += 1

	thead = """
<thead>
<tr>"""
	for sample in sorted(samples):
		thead = thead + "<th>{}</th>".format(sample)
		thead = thead + """
</tr>
</thead>"""
		for trans in sorted(transcript_counts):
			for sample in sorted(transcript_counts[trans]):
				thead = thead + """
<tbody>
<tr>
<td>{0}</td>
<td>{1}</td>
<td>{2}</td>
<td><a href="{3}">Link</a></td>
</tr>"""
	table = thead
	return table