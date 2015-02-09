#!/usr/bin/python

########################################################################
# 31 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

def ending():
	ending = """
		</tbody>
</table>
</table> 
</div>
</body>
</html>"""
	return ending

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
<li><a href="counts.html">Counts</a></li>
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
<li><a href="counts.html">Counts</a></li>
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
<div class="container" style="height: 600px; overflow: auto;">
<table class="table table-striped table-bordered table-condensed">
			"""
	return header

def p_header():
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
        <h1>Pynoncode</h1>
        <p><a href="#" class="btn btn-primary btn-lg" role="button">Learn more about pynoncode &raquo;</a></p>
      </div>
</div>
      """
	return header


def counts_trans_header():
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
<li <a href="pynoncode.html">Home</a></li>
<li class="active"><a href="#">Counts</a></li>
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