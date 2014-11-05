#!/usr/bin/python

########################################################################
# 31 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

def t_header():
	header="""
<!DOCTYPE html>
<html>
  <head>
    <title>Pynoncode</title>   
     <link href="dist/css/bootstrap.min.css" rel="stylesheet" media="screen">
    <link href="dist/css/bootstrap.min.css" rel="dist/css/bootstrap-theme.min.css">
    <script src="http://code.jquery.com/jquery-latest.js"></script>
    <script src="dist/js/bootstrap.min.js"></script>
    <meta name="viewport" content="width=device-width, initial-scale=1.0"> 
</head>
  <body>


<!-- Fixed navbar -->
<div class="navbar navbar-default" role="navigation">
<div class="container">
<div class="navbar-header">
<button type="button" class="navbar-toggle" data-toggle="collapse" data-target=".navbar-collapse">
<span class="sr-only">Toggle navigation</span>
<span class="icon-bar"></span>
<span class="icon-bar"></span>
<span class="icon-bar"></span>
        
</button>
<a class="navbar-brand" href="#">pynoncode</a>
</div>
<div class="navbar-collapse collapse">
<ul class="nav navbar-nav">
<li class="active"><a href="#">Home</a></li>
<li><a href="#about">About</a></li>
<li><a href="#contact">Contact</a></li>
</ul>
</div><!--/.nav-collapse -->
</div>
</div>

<div class="container"  role="main">

<div class="jumbotron">
        <h1>Pynoncode Transcripts results</h1>
        <p><a href="#" class="btn btn-primary btn-lg" role="button">Learn more about pynoncode &raquo;</a></p>
      </div>
<table class="table table-striped table-bordered table-condensed">
      
      """
	return header

def f_header():
	header="""
<!DOCTYPE html>
<html>
  <head>
    <title>Pynoncode</title>   
     <link href="dist/css/bootstrap.min.css" rel="stylesheet" media="screen">
    <link href="dist/css/bootstrap.min.css" rel="dist/css/bootstrap-theme.min.css">
    <script src="http://code.jquery.com/jquery-latest.js"></script>
    <script src="dist/js/bootstrap.min.js"></script>
    <meta name="viewport" content="width=device-width, initial-scale=1.0"> 
</head>
  <body>


<!-- Fixed navbar -->
<div class="navbar navbar-default" role="navigation">
<div class="container">
<div class="navbar-header">
<button type="button" class="navbar-toggle" data-toggle="collapse" data-target=".navbar-collapse">
<span class="sr-only">Toggle navigation</span>
<span class="icon-bar"></span>
<span class="icon-bar"></span>
<span class="icon-bar"></span>
        
</button>
<a class="navbar-brand" href="#">pynoncode</a>
</div>
<div class="navbar-collapse collapse">
<ul class="nav navbar-nav">
<li class="active"><a href="#">Home</a></li>
<li><a href="#about">About</a></li>
<li><a href="#contact">Contact</a></li>
</ul>
</div><!--/.nav-collapse -->
</div>
</div>

<div class="container"  role="main">

<div class="jumbotron">
        <h1>Pynoncode Fragments results</h1>
        <p><a href="#" class="btn btn-primary btn-lg" role="button">Learn more about pynoncode &raquo;</a></p>
      </div>
</div>
<div class="container" style="height: 800px; overflow: auto;">
<table class="table table-striped table-bordered table-condensed">
      """
	return header

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
	for trans in transcript_info:
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
		for trans in transcipts:
			for frag_pos in transcipts[trans]: #Now a list containing fragment, lfc, pvalue, padj
				link = "plots/{}.png".format(trans)
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
 	  				transcipts[trans][frag_pos][2][0], transcipts[trans][frag_pos][2][1], trans, link)
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
		for trans in transcipts:
			for frag_pos in transcipts[trans]: #Now a list containing fragment, lfc, pvalue, padj
				link = "plots/{}.png".format(trans)
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
 	  			""".format(frag_pos[0], frag_pos[1], frag_pos[2], transcipts[trans][frag_pos][0], transcipts[trans][frag_pos][1][0], transcipts[trans][frag_pos][1][1], trans, link)
	return thead

def ending():
	ending = """
		</tbody>
</table>
</table> 
</div>
</body>
</html>"""
	return ending

def create_transcript_html(transcripts):
	total_string = t_header()
	total_string += trans_table(transcripts)
	total_string += ending()
	return total_string

def create_fragment_html(transcripts, paired):
	total_string = f_header()
	total_string += frags_table(transcripts, paired)
	total_string += ending()
	return total_string