#!/bin/python
'''
Calculate the RNA-seq reads coverage over gene body.

Note:
1) Only input sorted and indexed BAM file(s). SAM format is not supported.
2) Genes/transcripts with mRNA length < 100 will be skipped (Number specified to "-l" cannot be < 100).
'''

#import built-in modules
import os,sys
if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import collections
import math
from time import strftime
import subprocess
from os.path import basename
import operator
#import third-party modules
from numpy import std,mean
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
import pysam

#import my own modules
from qcmodule import getBamFiles
from qcmodule import mystat
#changes to the paths
import pandas as pd
#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def valid_name(s):
	'''make sure the string 's' is valid name for R variable'''
	symbols = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'
	digit = '0123456789'
	rid = '_'.join(i for i in s.split())	#replace space(s) with '_'
	if rid[0] in digit:rid = 'V' + rid
	tmp = ''
	for i in rid:
		if i in symbols:
			tmp = tmp + i
		else:
			tmp = tmp + '_'
	return tmp


def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('log.txt','a')
	print(mesg, file=sys.stderr)
	print(mesg, file=LOG)


def pearson_moment_coefficient(lst):
	'''measure skewness'''
	mid_value = lst[int(len(lst)/2)]
	sigma = std(lst, ddof=1)
	tmp = []
	for i in lst:
		tmp.append(((i - mid_value)/sigma)**3)
	return mean(tmp)

def genebody_percentile(refbed, mRNA_len_cut = 100):
	'''
	return percentile points of gene body
	mRNA length < mRNA_len_cut will be skipped
	'''
	if refbed is None:
		print("You must specify a bed file representing gene model\n", file=sys.stderr)
		exit(0)

	g_percentiles = {}
	transcript_count = 0
	for line in open(refbed,'r'):
		try:
			if line.startswith(('#','track','browser')):continue
			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5]
			geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])

			exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
			exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
			exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
			exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
			transcript_count += 1
		except:
			print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
			continue
		gene_all_base=[]
		mRNA_len =0
		flag=0
		for st,end in zip(exon_starts,exon_ends):
			gene_all_base.extend(list(range(st+1,end+1)))		#1-based coordinates on genome
		if len(gene_all_base) < mRNA_len_cut:
			continue
		g_percentiles[geneID] = (chrom, strand, mystat.percentile_list (gene_all_base))	#get 100 points from each gene's coordinates
	printlog("Total " + str(transcript_count) + ' transcripts loaded')
	return g_percentiles

def genebody_coverage(bam, position_list):
	'''
	position_list is dict returned from genebody_percentile
	position is 1-based genome coordinate
	'''
	samfile = pysam.Samfile(bam, "rb")
	aggreagated_cvg = collections.defaultdict(int)
	gene_finished = 0
	# by Raza
	string = "percent_"
	cols=[string+str(i) for i in list(range(1, 101))]
	cols.insert(0,"GeneName")
	lst = []
	for key,value in position_list.items():
		chrom=value[0]
		strand=value[1]
		positions=value[2]
		coverage = {}
		for i in positions:
			coverage[i] = 0.0
		chrom_start = positions[0]-1
		if chrom_start <0: chrom_start=0
		chrom_end = positions[-1]
		try:
			samfile.pileup(chrom, 1,2)
		except:
			continue

		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in positions:
				continue
			if pileupcolumn.n == 0:
				coverage[ref_pos] = 0
				continue
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue
				if pileupread.alignment.is_secondary:continue
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read +=1
			coverage[ref_pos] = cover_read
		tmp = [coverage[k] for k in sorted(coverage)]
		if strand == '-':
			tmp = tmp[::-1]
		for i in range(0,len(tmp)):
			aggreagated_cvg[i] += tmp[i]

		if len(tmp)==99 :
			tmp.append(0)

		tmp.insert(0,key.split('_')[3])
		lst.append(tmp)
		gene_finished += 1

		if gene_finished % 100 == 0:
			print("\t%d transcripts finished\r" % (gene_finished), end=' ', file=sys.stderr)
	df = pd.DataFrame(lst, columns=cols)
	return 	aggreagated_cvg,df

def Rcode_write(dataset,file_prefix, format='pdf', colNum=100):
	'''generate R script for visualization'''
	ROUT = open(file_prefix + '.r','w')
	names=[]
	datas=[]
	for name, data, tmp in dataset:
		names.append(name)
		datas.append(data)
		print(name + ' <- c(' + ','.join([str(i) for i in data]) + ')', file=ROUT)

	tick_pos = [1,10,20,30,40,50,60,70,80,90,100]
	tick_lab = [1,10,20,30,40,50,60,70,80,90,100]

	# do not generate heatmap if only 1 sample
	if len(names) >=3:
		print('data_matrix' + ' <- matrix(c(' + ','.join(names) + '), byrow=T, ' +  'ncol=' + str(colNum) + ')', file=ROUT)
		print('rowLabel <- c(' + ','.join(['"' + i + '"' for i in names]) + ')', file=ROUT)
		print('\n', file=ROUT)
		print('%s(\"%s.%s\")' % (format.lower(),file_prefix + ".heatMap",format.lower()), file=ROUT)
		print('rc <- cm.colors(ncol(data_matrix))', file=ROUT)
		print('heatmap(data_matrix' + ', scale=c(\"none\"),keep.dendro=F, labRow = rowLabel ' + ',Colv = NA,Rowv = NA,labCol=NA,col=cm.colors(256),margins = c(6, 8),ColSideColors = rc,cexRow=1,cexCol=1,xlab="Gene body percentile (5\'->3\')", add.expr=x_axis_expr <- axis(side=1,at=c(%s),labels=c(%s)))' % (','.join([str(i) for i in tick_pos]), ','.join(['"' + str(i) + '"' for i in tick_lab])), file=ROUT)
		print('dev.off()', file=ROUT)


	print('\n', file=ROUT)

	print('%s(\"%s.%s\")' % (format.lower(),file_prefix + ".curves",format.lower()), file=ROUT)
	print("x=1:%d" % (colNum), file=ROUT)
	print('icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(%d)' % (len(names)), file=ROUT)

	if len(names) == 1:
		print("plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0]), file=ROUT)

	elif  len(names) >=2 and len(names) <=6:
		print("plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0]), file=ROUT)
		for i in range(1,len(names)):
			print("lines(x,%s,type='l',col=icolor[%d])" % (names[i], i+1), file=ROUT)
		print("legend(0,1,fill=icolor[%d:%d], legend=c(%s))" % (1,len(names), ','.join([ "'" + str(n) + "'" for n in names])), file=ROUT)

	elif len(names) > 6:
		print('layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), 4, 4, byrow = TRUE))', file=ROUT)
		print("plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0]), file=ROUT)
		for i in range(1,len(names)):
			print("lines(x,%s,type='l',col=icolor[%d])" % (names[i], i+1), file=ROUT)
		print('par(mar=c(1,0,2,1))', file=ROUT)
		print('plot.new()', file=ROUT)
		print("legend(0,1,fill=icolor[%d:%d],legend=c(%s))" % (1,len(names), ','.join([ "'" + str(n) + "'" for n in names])), file=ROUT)

	print('dev.off()', file=ROUT)
	ROUT.close()
def writePerGeneCoverage(df,file_prefix):
	df.to_csv(file_prefix, encoding='utf-8',index=False)

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_files",help='Input file(s) in BAM format. "-i" takes these input: 1) a single BAM file. 2) "," separated BAM files. 3) directory containing one or more bam files. 4) plain text file containing the path of one or more bam file (Each row is a BAM file path). All BAM files should be sorted and indexed using samtools.')
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. [required]")
	parser.add_option("-l","--minimum_length",action="store",type="int",default=100, dest="min_mRNA_length",help="Minimum mRNA length (bp). mRNA smaller than \"min_mRNA_length\" will be skipped. default=%default")
	parser.add_option("-f","--format",action="store",type="string",dest="output_format", default='pdf', help="Output file format, 'pdf', 'png' or 'jpeg'. default=%default")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_files and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)

	if not os.path.exists(options.ref_gene_model):
		print('\n\n' + options.ref_gene_model + " does NOT exists" + '\n', file=sys.stderr)
		#parser.print_help()
		sys.exit(0)
	if options.min_mRNA_length < 100:
		print('The number specified to "-l" cannot be smaller than 100.' + '\n', file=sys.stderr)
		sys.exit(0)

	OUT1 = open(options.output_prefix  + ".geneBodyCoverage.txt"	,'w')
	print("Percentile\t" + '\t'.join([str(i) for i in range(1,101)]), file=OUT1)

	printlog("Read BED file (reference gene model) ...")
	gene_percentiles = genebody_percentile(refbed = options.ref_gene_model, mRNA_len_cut = options.min_mRNA_length)

	printlog("Get BAM file(s) ...")
	bamfiles = getBamFiles.get_bam_files(options.input_files)
	for f in bamfiles:
		print("\t" + f, file=sys.stderr)

	file_container = []
	for bamfile in bamfiles:
		printlog("Processing " + basename(bamfile) + ' ...')
		cvg,df = genebody_coverage(bamfile, gene_percentiles)
		if len(cvg) == 0:
			print("\nCannot get coverage signal from " + basename(bamfile) + ' ! Skip', file=sys.stderr)
			continue
		tmp = valid_name(basename(bamfile).replace('.bam',''))	# scrutinize R identifer
		writePerGeneCoverage(df, options.output_prefix + '.perGeneBodyCoverage.txt')
		if file_container.count(tmp) == 0:
			print(tmp + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)]), file=OUT1)
		else:
			print(tmp + '.' + str(file_container.count(tmp)) + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)]), file=OUT1)
		file_container.append(tmp)
	OUT1.close()


	dataset=[]
	for line in open(options.output_prefix  + ".geneBodyCoverage.txt",'r'):
		line = line.strip()
		if line.startswith("Percentile"):
			continue
		f = line.split()
		name = f[0]
		dat = [float(i) for i in  f[1:]]
		skewness = pearson_moment_coefficient(dat)
		dataset.append((name, [(i -min(dat))/(max(dat) - min(dat)) for i in dat], skewness))
	dataset.sort(key = operator.itemgetter(2), reverse=True)

	print("\n\n", file=sys.stderr)
	print("\tSample\tSkewness", file=sys.stderr)
	for a,b,c in dataset:
		print('\t' + a + '\t' + str(c), file=sys.stderr)
	Rcode_write(dataset, options.output_prefix + '.geneBodyCoverage', format = options.output_format)

	printlog("Running R script ...")
	try:
		subprocess.call("Rscript " + options.output_prefix + '.geneBodyCoverage.r',shell=True)
	except:
		print("Cannot generate pdf file from " + options.output_prefix + '.geneBodyCoverage.r', file=sys.stderr)
		pass


if __name__ == '__main__':
	main()
