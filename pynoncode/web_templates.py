#!/usr/bin/python

########################################################################
# 31 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

def header():
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

def ending():
	ending = """
		</tbody>
</table>
</div> 
</body>
</html>"""
	return ending

def create_html(transcripts, ):
	total_string = header()
	total_string += trans_table(transcripts)
	total_string += ending()
	return total_string