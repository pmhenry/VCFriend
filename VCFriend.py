#!/usr/bin/env python3
#
# VCFriend
#
######################################################################
#
# VCFriend is a collection of functions for manipulating and
# analyzing data from VCF files. 
#
######################################################################
import argparse
from apps import utilities
from apps import allele_seq, clean, compare, matrix
from apps import pat_match, sim_matrix, remove

#####################################################################################
# Application Subparsers

###############
# allele_seq
  
def allele_seqParser(subparsers):
	allele_seq_parser = subparsers.add_parser('allele-seq',
		help='Creates a multi-fasta file containing the sequence of all snp alleles for each sample in a VCF (haploid and snps only for now)')
	allele_seq_parser.add_argument('-i', '--input', help='VCF file', dest='InFile', type=str, default='Error1')
	allele_seq_parser.add_argument('-o', '--output', help='Name of output file', dest='OutFile', type=str, default='Error2')

	return allele_seq_parser

class allele_seqCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = allele_seq.allele_seqApp()
   		return app.start(args.InFile, args.OutFile)

###############
# clean
 
def cleanParser(subparsers):
	clean_parser = subparsers.add_parser('clean',
		help='Removes all variants from VCF files or a Variant Matrix where a genotype call could not be made for any isolate')
	clean_parser.add_argument('-i', '--input', help='VCF or Table file', dest='InFile', type=str, default='Error1')
	clean_parser.add_argument('-o', '--output', help='Name of output table', dest='OutFile', type=str, default='Error2')
	
	return clean_parser

class cleanCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = clean.cleanApp()
   		return app.start(args.InFile, args.OutFile)


 ###############
# compare
  
def compareParser(subparsers):
	compare_parser = subparsers.add_parser('compare',
		help='Finds the percent of variants shared or not between samples in a VCF file or Variant Matrix')
	compare_parser.add_argument('-i', '--input', help='VCF Table (output from vcf_mat)', dest='InFile', type=str, default='Error1')
	compare_parser.add_argument('-s', '--samples', help='List of samples to include separated by comma', dest='Samples', type=str, default='Error2')
	compare_parser.add_argument('-x', '--exclude', help='List of samples to exclude separated by comma (optional)', dest='Exclude', type=str, default=None)
	compare_parser.add_argument('-o', '--output', help='output file for tagged variants (optional)', dest='OutFile', type=str, default=None)
	
	return compare_parser

class compareCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = compare.compareApp()
   		return app.start(args.InFile, args.Samples, args.Exclude, args.OutFile)


###############
# matrix 

def matrixParser(subparsers):
	matrix_parser = subparsers.add_parser('matrix', 
		help='Converts VCF file into tab deliminated text file where samples are rows and variants are columns (VCF Matrix)')
	matrix_parser.add_argument('-i', '--intput', help='VCF file', dest='InFile', type=str, default='Error1')
	matrix_parser.add_argument('-o', '--output', help='Name of output table', dest='OutFile', type=str, default='Error2')
	
	return matrix_parser

class matrixCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = matrix.matrixApp()
   		return app.start(args.InFile, args.OutFile)


###############
# pat_match 

def pat_matchParser(subparsers):
	pat_match_parser = subparsers.add_parser('pat-match', 
		help='Extracts variants matching designated patterns of presence and absence from VCF files or Variant Matrix')
	pat_match_parser.add_argument('-i', '--input', help='VCF File or Table (output from vcf_mat)', dest='InFile', type=str, 
		default='Error1')
	pat_match_parser.add_argument('-p', '--pattern', help='Comma separated string in quotes containing genotype calls or N\'s, Y\'s, and .\'s. Pattern is phase sensitive.',
		dest='Pattern', type=str, default='Error2')
	pat_match_parser.add_argument('-o', '--output', help='Name of output file', dest='OutFile', type=str, default='Error3')
	
	return pat_match_parser

class pat_matchCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = pat_match.pat_matchApp()
   		return app.start(args.InFile, args.Pattern, args.OutFile)

###############
# remove

def removeParser(subparsers):
	remove_parser = subparsers.add_parser('remove', 
		help=' Removes samples from a VCF or Variant Matrix file.')
	remove_parser.add_argument('-i', '--input', help='VCF or Table file', dest='InFile', type=str, default='Error1')
	remove_parser.add_argument('-s', '--sample', help='Comma separated list of samples in quotes', dest='Sample', type=str, default='Error2')
	remove_parser.add_argument('-o', '--output', help='Name of output VCF', dest='OutFile', type=str, default='Error3')
	
	return remove_parser

class removeCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = remove.removeApp()
   		return app.start(args.InFile, args.Sample, args.OutFile)


###############
# sim-matrix
  
def sim_matrixParser(subparsers):
	sim_matrix_parser = subparsers.add_parser('sim-matrix',
		help='Create a similarity matrix of all samples in a a VCF file or Variant Matrix')
	sim_matrix_parser.add_argument('-i', '--input', help='VCF or Table File (output from matrix)', dest='InFile', type=str, default='Error1')
	sim_matrix_parser.add_argument('-o', '--output', help='Name of output similarity matrix', dest='OutFile', type=str, default='Error2')

	return sim_matrix_parser

class sim_matrixCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = sim_matrix.sim_matrixApp()
   		return app.start(args.InFile, args.OutFile)


#####################################################################################
# Program Parser

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Functions for manipulating and analyzing VCF files.', add_help=True,
        epilog="For questions or comments, contact Bradley Jenner <bnjenner@ucdavis.edu>")
    subparsers = parser.add_subparsers(help='commands', dest='command')
    allele_seqParser(subparsers)
    cleanParser(subparsers)
    compareParser(subparsers)
    matrixParser(subparsers)
    pat_matchParser(subparsers)
    removeParser(subparsers)
    sim_matrixParser(subparsers)
    
    args = parser.parse_args()
    return args

def main():
	allele_seq = allele_seqCMD()
	clean = cleanCMD()
	compare = compareCMD()
	matrix = matrixCMD()
	pat_match = pat_matchCMD()
	remove = removeCMD()
	sim_matrix = sim_matrixCMD()
	
	commands = {'allele-seq': allele_seq, 'clean': clean, 'compare': compare, 'matrix': matrix, 
				'pat-match': pat_match, 'remove': remove, 'sim-matrix': sim_matrix 
			   } 

	args = parseArgs()
	commands[args.command].execute(args)
	
if __name__ == '__main__':
	main()
