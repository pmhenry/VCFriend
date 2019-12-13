# VCFriend
#
######################################################################
#
# VCFriend is a collection of functions for manipulating and
# analyzing data from VCF files. 
#
######################################################################
import os
import argparse

class removeApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, Sample, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF File Required*')
		elif Sample == 'Error2':
			raise Exception('*Sample ID Required*')
		elif OutFile == 'Error3':
			raise Exception('*Name for Output File Required*')

		hash_lines = []
		var_lines =[]

		with open(InFile, 'r') as fi:  # read in VCF file 
			vcf = fi.readlines()
			for line in vcf:
				if '#' in line:
					hash_lines.append(line)  # appends all header lines to hash_lines list
				else:
					var_lines.append(line)  # appends all non-header lines to var_lines list


		header = hash_lines[-1]  # takes header lines in vcf file and removes it from original list.
		hash_lines.remove(hash_lines[-1])

		for x in range(len(header.split('\t'))):  # iterates over list where each item is a column 
			if header.split('\t')[x] == Sample:   # takes column number that matches sample ID to remove
				rm_1 = x

		final_vars = []  
		temp_head = header.split('\t')
		temp_head.remove(header.split('\t')[rm_1])  # removes isolate from header line 
		final_head = '\t'.join(temp_head)


		for line in var_lines:  # iterates over every line in var_lines list 
			temp = line.split('\t')  # creates temp list where each item is a column from the VCF file
			temp.remove(temp[rm_1])  # removes column corresponding to removed sample
			temp_2 = '\t'.join(temp)
			final_vars.append(temp_2)


		with open(OutFile, 'w') as fo:  # creates file with original hash lines from input VCF.
			for line in hash_lines:
				fo.write(line)
			fo.write(final_head)  # addition of new headers
			for line in final_vars:
				fo.write(line)  # addition of new variant lines


class vcf_matApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Name for Output File Required*')

		vcf = []
		snp_mat = []

		with open(InFile, 'r') as fi:  # read in VCF file and appends all non-header lines 
			file = fi.readlines()      # to a list-- 'vcf'. Each item is a new line.
			for line in file:
				if '##' not in line:
					vcf.append(line)
		
		isolates = vcf[0].split('\t')[9:] # extracts isolate names from VCF. List contains isolate names in order
		isolates[-1] = isolates[-1][:-1] # removes newline at end of string
		vcf.pop(0) # removes header from VCF data
		head = ['']

		for x in range(len(isolates)): # cycles through isolate columns list
			temp = [isolates[x]]

			for y in range(len(vcf)):  # cycles through each line in VCF file.
				temp.append(vcf[y].split('\t')[9+x].split(':')[0].strip())  # takes genotype assignment for isolate on each line
			
				if len(head) <= len(vcf):
					name = str(vcf[y].split('\t')[0])+':'+str(vcf[y].split('\t')[1]) # extracts variants scaffold and position
					head.append(name) # unique variant name is 'scaffold':'position'
			
			snps = '\t'.join(temp)
			snp_mat.append(snps)  # appends string containing genotype calls separated by tabs  
		

		header = '\t'.join(head) # creates tab separated string containing variant ID's.

		with open(OutFile, 'w') as fo:  # creates tab separated tables where columns are variants and rows are isolate genotypes
			fo.write(header+'\n')
			for row in snp_mat:
				fo.write(row+'\n')


class pat_matchApp():

	def __init__(self):
		self.verbose = False

	def pat_match(pattern, text):  #  matching algorithm, returns either '1' (match) or '0' (no match)  

			for i in range(len(pattern)):  #  compares pattern and text by individual characters
				match = True
				if pattern[i] == text[i] or pattern[i] == 'N':  #  loop continues and match remains equal to true if pattern matches  
					continue                                    #  text at given position or pattern at the given position is 'N'
				elif pattern[i] == 'Y' and text[i] != '.':
					continue                               
				else:  #  if pattern does not equal match at given position, match is set to false and breaks out of the loop
					match = False
					break		
			if match == True:
				m = 1
			else:
				m = 0

			return m

	def start(self, InFile, Pattern, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF Table File Required*')
		elif Pattern == 'Error2':
			raise Exception('*Presence/Absence Pattern Required*')
		elif OutFile == 'Error3':
			raise Exception('*Name for Output File Required*')


		with open(InFile, 'r') as fi:
			file = fi.readlines()
			fi.close()

		l = len(file[0].split('\t'))  # length variable for loop below
		gene_lst = file[0].split('\t')[1:]  # list of gene headers in order given by file
		

		text_lst = []

		for x in range(1, l):  # parses tab deliminated text file for pattern and generates a list of texts to search 
			temp = []
			for line in file:  # splits ever line of file into a list by tabs, assigns each item on each line at a given position to a list
				number = line.split('\t')[x] 
				if number == '':
					temp.append('.')
				elif number[-1] == '\n':  # removes \n characters
					temp.append(number[:-1])
				else:
					temp.append(number)
			temp.pop(0)  # removes title of each column 
			text_lst.append(''.join(temp)) 

		match_lst = []
		matches = []

		for text in text_lst:  # loop for matching function that compares pattern to each text in the list of texts
			match_lst.append(pat_matchApp.pat_match(Pattern, text))  # outputs a '1' or '0' and appends to a list 

		# the final result of the above loop is a list of '0's and '1's whose position in the list 
		# matches the the position of its corresponding gene in the list of genes	

		for x in range(len(gene_lst)):  # indexes list of matches for instances of '1' and extracts the corresponding gene header at the
		    if match_lst[x] == 1:       # same position in the list of gene headers.  Assigns genes to final lisst (matches) 
		    	matches.append(gene_lst[x])


		with open(OutFile, 'w') as fo:  # creates a new file with every line being a gene header + '\n'
		    for gene in matches:
		    	fo.write(gene.strip()+'\n')
		    fo.close()


class nono_callsApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF Table File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Name for Output File Required*')

		hash_lines = []
		var_lines =[]

		with open(InFile, 'r') as fi: # reads in VCF table file 
			vcf = fi.readlines()  
			for line in vcf:
				if '#' in line:
					hash_lines.append(line)  # appends all header lines to hash_lines list
				else:
					var_lines.append(line)  # appends all non-header lines to var_lines list
		

		no_calls = []

		for x in range(len(var_lines)):  # iterates over all lines in var_lines list  
			for y in range(9, len(var_lines[0].split('\t'))):  # iterates over every volumn in every variant line
				if '.' in var_lines[x].split('\t')[y][0]: 
					no_calls.append(x)  # adds line # to no_calls list if no genotype call is detected
					break

		no_calls = list(dict.fromkeys(no_calls))  # reverses order of lines in list 
		no_calls.sort(reverse=True)


		for site in no_calls:  # iterates over no_calls list and deletes the corresponding line in var_lines list
			del var_lines[site]

		with open(OutFile, 'w') as fo:  # creates file with same hash lines as input file and only variant lines without no calls
			for line in hash_lines:
				fo.write(line)
			for line in var_lines:
				fo.write(line) 


class samp_compApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, Samples, Exclude, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF Table File Required*')
		elif Samples == 'Error2':
			raise Exception('*Samples Required*')

		with open(InFile, 'r') as fi: # reads in VCF table from vat_mat command 
			file = fi.readlines()

			
		if len(Samples.split(',')) < 2:
			raise Exception('*More Than One Sample Required*')
		else:
			samples = Samples.split(',')

		ex_list = []

		try:	
			ex_samples = Exclude.split(',')
			print('Excluded Samples: ', end="")
			print(ex_samples)
		except:
			ex_samples = []

		
		print('Total Variants: ' + str(int(len(file[0].split('\t')))-1))  # counts number of columns 

		samples = set(Samples.split(','))

		if len(samples) == 1 and len(ex_samples) == 0:
			print('Shared Variants: ' + str(int(len(file[0].split('\t')))-1))
			print('Similarity: 1.0') 
			return

		samp_list = []

		for i in range(len(file)):  # iterates over every row
			if file[i].split('\t')[0] in samples:  # finds line with corresponding first sample ID
				samp_list.append(file[i].split('\t')[1:])  # creates list with every genotype call for every variant for first isolate
				if len(samples) == 1:
					samp_list.append(file[i].split('\t')[1:])

			elif file[i].split('\t')[0] in ex_samples:
				ex_list.append(file[i].split('\t')[1:])
		

		sim = 0 
		sim_list = []

		for x in range(len(samp_list[0])):
			add = True
			for y in range(1, len(samp_list)):
				if samp_list[0][x] == samp_list[y][x]:
					pass
				else:
					add = False
					break
			if add == True:
				for z in range(len(ex_list)):
					if samp_list[0][x] == ex_list[z][x]:
						add = False
						break
				if add == True:
					sim += 1
					sim_list.append(file[0].split('\t')[x+1]) 


		if Exclude == None:
			print('Shared Variants: ' + str(sim)) 
		else:
			print('Variants not shared with excluded samples: ' + str(sim))


		if OutFile != None:
			with open(OutFile, "w") as fo:
				for sim in sim_list:
					fo.write(sim + '\n')


class dist_matApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF Table File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Output File Name Required*')

		table = {}
		samples = []
		dist = []

		with open(InFile, "r") as fi:  # reads in vcf matrix
			lines = fi.readlines()

		for line in lines:  # creates dictionary with sample as key an allele sequence as value 
			table[line.split("\t")[0]] = "".join(line.split("\t")[1:])
			samples.append(line.split("\t")[0])  # also adds all samples to a list


		num_samps = len(samples) 

		for i in range(num_samps - 1):  # populates the distance matrix with 0's to the propper dimensions
			dist.append([])
			for j in range(num_samps - 1):
				dist[i].append("0")

		for x in range(1, num_samps):  # iterates over columns
			for y in range(1, num_samps): # iterates over columns for comparison
				sim = 0
				samp_1 = table[samples[x]]
				samp_2 = table[samples[y]]
				if len(samp_1) >= len(samp_2):  # choses shortest allele sequence
					length = len(samp_1)
				else:
					length = len(samp_2)
				for c in range(length):
					if samp_1[c] == samp_2[c]:  # for each match in allele sequences, similarity count increases
						sim += 1
				perc_sim = sim / (( len(samp_1) + len(samp_2) ) / 2 )  # percent similarity formula
				dist[x-1][y-1] = perc_sim  # adds above value to correct spot in matrix

		with open(OutFile, "w") as fo:  # opens output file
			fo.write('\t'.join(samples) + '\n')  # populates matrix column headers
			for x in range(len(dist)):           # populates matrix by row
				temp = samples[x+1] 
				for y in range(len(dist[x])):
					temp += '\t' + str(dist[x][y])
				fo.write(temp + '\n') 



class snp_statApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Output File Suffix Required*')

		input_file = InFile
		output_prefix = '/'.join(OutFile.split('/')[:-1]) + '/'
		output_suffix = OutFile.split('/')[-1]
		var_output = output_prefix + "variant_" + output_suffix
		scaff_output = output_prefix + "scaffold_" + output_suffix

		header = []
		variant_list = []

		scaff_dict = {}
		metrics_dict = {
						"DP" : 0.0,
						"QD" : 0.0,
						"SOR" : 0.0,
						"MQ" : 0.0,
						"MAPQ" : 0.0,
						"MQRankSum" : 0.0,
						"BaseQRankSum" : 0.0
						}

		metrics_keys = metrics_dict.keys()


		with open(input_file, "r") as fi:
			lines = fi.readlines()

			for line in lines:

				if "#CHROM" in line:
					header = line.split('\t')[:7]
				elif "#" not in line:
					variant_list.append(line.split('\t'))

			header[0] = "CHROM"


		with open(var_output, "w") as fo_v:

			fo_v.write('\t'.join(header) + '\n')


			for i in range(len(variant_list)):

				temp_list = variant_list[i][:7]
				var_info_field = variant_list[i][7]

				scaffold = variant_list[i][0]

				if scaffold not in scaff_dict.keys():
					scaff_dict[scaffold] = [1, metrics_dict.copy()]

				else: 
					scaff_dict[scaffold][0] += 1

				for key in metrics_keys:

					for info in var_info_field.split(';'):

						if info.split('=')[0] == key and info.split('=')[0] != info.split('=')[-1]:

							scaff_dict[scaffold][1][key] += float(info.split('=')[-1])
							temp_list.append(info.split('=')[-1])
							break
							

					else:
						temp_list.append("NA")


			
				fo_v.write('\t'.join(temp_list) + '\n')


		with open(scaff_output, "w") as fo_s:

			fo_s.write('CHROM' + '\t' + '\t'.join(metrics_keys) + '\n')

			for scaff_key in scaff_dict.keys():
				for met_key in metrics_keys:
					scaff_dict[scaff_key][1][met_key] = str(round((float(scaff_dict[scaff_key][1][met_key]) / scaff_dict[scaff_key][0]), 3))

				fo_s.write(scaff_key + '\t' + '\t'.join(scaff_dict[scaff_key][1].values()) + '\n')

###############
# remove

def removeParser(subparsers):
	remove_parser = subparsers.add_parser('remove', 
		help='Removes sample and its corresponding column from VCF')
	remove_parser.add_argument('-i', '--input', help='VCF file', dest='InFile', type=str, default='Error1')
	remove_parser.add_argument('-s', '--sample', help='Sample ID for Removal', dest='Sample', type=str, default='Error2')
	remove_parser.add_argument('-o', '--output', help='Name of output VCF', dest='OutFile', type=str, default='Error3')
	
	return remove_parser

class removeCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = removeApp()
   		return app.start(args.InFile, args.Sample, args.OutFile)

###############
# vcf_mat 

def vcf_matParser(subparsers):
	vcf_mat_parser = subparsers.add_parser('vcf_mat', 
		help='Converts VCF file into tab deliminated text file where samples are rows and variants are columns')
	vcf_mat_parser.add_argument('-i', '--intput', help='VCF file', dest='InFile', type=str, default='Error1')
	vcf_mat_parser.add_argument('-o', '--output', help='Name of output table', dest='OutFile', type=str, default='Error2')
	
	return vcf_mat_parser

class vcf_matCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = vcf_matApp()
   		return app.start(args.InFile, args.OutFile)

###############
# pat_match 

def pat_matchParser(subparsers):
	pat_match_parser = subparsers.add_parser('pat_match', 
		help='Extracts variants matching designated patterns of presence and absence.')
	pat_match_parser.add_argument('-i', '--input', help='VCF Table (output from vcf_mat)', dest='InFile', type=str, 
		default='Error1')
	pat_match_parser.add_argument('-p', '--pattern', help='String containing only 1\'s, 0\'s, N\'s, Y\'s, and .\'s.',
		dest='Pattern', type=str, default='Error2')
	pat_match_parser.add_argument('-o', '--output', help='Name of output file', dest='OutFile', type=str, default='Error3')
	
	return pat_match_parser

class pat_matchCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = pat_matchApp()
   		return app.start(args.InFile, args.Pattern, args.OutFile)

###############
# nono_calls
 
def nono_callsParser(subparsers):
	nono_calls_parser = subparsers.add_parser('nono_calls',
		help='Removes all variants where a genotype call could not be made for any isolate')
	nono_calls_parser.add_argument('-i', '--input', help='VCF file', dest='InFile', type=str, default='Error1')
	nono_calls_parser.add_argument('-o', '--output', help='Name of output table', dest='OutFile', type=str, default='Error2')
	
	return nono_calls_parser

class nono_callsCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = nono_callsApp()
   		return app.start(args.InFile, args.OutFile)

###############
# samp_comp
  
def samp_compParser(subparsers):
	samp_comp_parser = subparsers.add_parser('samp_comp',
		help='Finds the percent of variants shared between two samples in a vcf table')
	samp_comp_parser.add_argument('-i', '--input', help='VCF Table (output from vcf_mat)', dest='InFile', type=str, default='Error1')
	samp_comp_parser.add_argument('-s', '--samples', help='List of samples to include separated by comma', dest='Samples', type=str, default='Error2')
	samp_comp_parser.add_argument('-x', '--exclude', help='List of samples to exclude separated by comma', dest='Exclude', type=str, default=None)
	samp_comp_parser.add_argument('-o', '--output', help='optional output file for tagged variants', dest='OutFile', type=str, default=None)
	return samp_comp_parser

class samp_compCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = samp_compApp()
   		return app.start(args.InFile, args.Samples, args.Exclude, args.OutFile)

###############
# dist_mat
  
def dist_matParser(subparsers):
	dist_mat_parser = subparsers.add_parser('dist_mat',
		help='Uses vcf table to create a distance matrix of all samples. NOTE: This is not optimized (every spot in distance matrix is calculated, even repeats)')
	dist_mat_parser.add_argument('-i', '--input', help='VCF Table (output from vcf_mat)', dest='InFile', type=str, default='Error1')
	dist_mat_parser.add_argument('-o', '--output', help='Name of output distance matrix', dest='OutFile', type=str, default='Error2')

	return dist_mat_parser

class dist_matCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = dist_matApp()
   		return app.start(args.InFile, args.OutFile)

###############
# snp_stat
  
def snp_statParser(subparsers):
	snp_statparser = subparsers.add_parser('snp_stat',
		help='Uses vcf file to create a brief summary table of certain stats at the variant and chromosome level. (Comments coming soon...)')
	snp_statparser.add_argument('-i', '--input', help='VCF File', dest='InFile', type=str, default='Error1')
	snp_statparser.add_argument('-o', '--output', help='Name of output suffix for summary tables', dest='OutFile', type=str, default='Error2')

	return snp_statparser

class snp_statCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = snp_statApp()
   		return app.start(args.InFile, args.OutFile)

#####################################################################################


def parseArgs():
    parser = argparse.ArgumentParser(
        description='Miscellaneous functions for manipulating and analyzing VCF files.', add_help=True,
        epilog="For questions or comments, contact Bradley Jenner <bnjenner@ucdavis.edu>")
    subparsers = parser.add_subparsers(help='commands', dest='command')
    removeParser(subparsers)
    vcf_matParser(subparsers)
    pat_matchParser(subparsers)
    nono_callsParser(subparsers)
    samp_compParser(subparsers)
    dist_matParser(subparsers)
    snp_statParser(subparsers)
    args = parser.parse_args()
    return args

def main():
	remove = removeCMD()
	vcf_mat = vcf_matCMD()
	pat_match = pat_matchCMD()
	nono_calls = nono_callsCMD()
	samp_comp = samp_compCMD()
	dist_mat = dist_matCMD()
	snp_stat = snp_statCMD()
	commands = {'remove':remove, 'vcf_mat':vcf_mat, 
				'pat_match':pat_match, 'nono_calls':nono_calls, 
				'samp_comp':samp_comp , 'dist_mat':dist_mat,
				'snp_stat':snp_stat}
	args = parseArgs()
	commands[args.command].execute(args)
	
if __name__ == '__main__':
	main()
