multiple scripts to be used in bioinformatics

Description
===========

base_query.py: used to investigate base positions in a BAM file

bb_fit.py: contains a python class for fitting a beta-binomial distribution

mixture_model.py: used to cluster allele frequencies using a gaussian mixture model

Invocation
==========

base_query.py::

	$ python base_query.py -h
	
Gives::

	optional arguments:
	  -h, --help     show this help message and exit
	  --bam BAM      input bam file
	  --chrom CHROM  query chromosome
	  --start START  query start
	  --end END      query end
	  --filter       flag to filter all regions
	  --summary      only print summary stuff


mixture_model.py::

	$ python mixture_model.py -h
	
Gives::

	optional arguments:
	  -h, --help  show this help message and exit
	  --vcf VCF   vcf file
	  --max MAX   max number of sub-clones


bb_fit.py::

	>>> from bb_fit import bb_model
	>>> observed = [1,2,3,2]
	>>> total = [10,9,11,12]
	>>> x = bb_model(observed, total) #fits beta-binomial, can take numpy arrays
	
	>>> x.mean #sample size 1
	0.19043834888300976
	>>> x.update_mean(10) #sample size 10
	>>> x.mean
	1.9043834888300974
	
	>>> x.var #sample size 1
	0.15417158415772284
	>>> x.update_var(10) #sample size 10
	>>> x.var
	0.6282503703278377
	
	>>> x.skew #sample size 1
	-1.5767940867704575
	>>> x.update_skew(10) #sample size 10
	>>> x.skew
	0.20973400768832143
	
	>>> x.alpha
	-3.0831764512182462
	
	>>> x.beta
	-13.106716337194376
	
	>>> x.cdf(1,10) #10 trials, 1 success
	0.2941167290438681
	
	>>> x.pdf(1,10) #10 trials, 1 success
	0.25954585054993118
	
	>>> x.plot_pdf(100) #sample size 100
	>>> x.plot_pdf(100, "fig1.png") #save as fig1.png
	
	>>> x.plot_cdf(100)#sample size 100
	>>> x.plot_cdf(100, "fig1.png") #save as fig1.png
	
OR::

	>>> x = bb_model(alpha=134, beta=145) #set alpha and beta
	>>> x.mean
	0.48028673835125446
	