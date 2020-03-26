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


###############
# remove

class removeApp():

	def __init__(self):
		self.verbose = False


	def start(self, InFile, Sample, OutFile):
	
		if InFile == 'Error1':
			raise Exception('*VCF or Table File Required*')
		elif Sample == 'Error2':
			raise Exception('*Sample ID Required*')
		elif OutFile == 'Error3':
			raise Exception('*Name for Output File Required*')

		var_lines =[]

		samp_list = Sample.split(',')

		with open(InFile, 'r') as fi:  # read in VCF file 
			file = fi.readlines()


		with open(OutFile, 'w') as fo:  # creates file with original hash lines from input VCF.

			if '##' in file[1]: # Checks if File is in VCF or Table format

				for line in vcf:
					if '##' in line:
						fo.write(line)
					else:
						if '#CHROM' in line:
							header = line.split('\t')
							for x in range(9, len(header)):
								if header[x] in samp_list:
									rm_1 = x # assigns index to remove in future lines
									header.remove(header[rm_1]) # removes column corresponding to removed sample
									fo.write('\t'.join(header))
									break
						else:
							temp = line.split('\t')  # creates temp list where each item is a column from the VCF file
							temp.remove(temp[rm_1])  # removes column corresponding to removed sample
							fo.write('\t'.join(temp))

			else:

				for line in file:
					clear = True
					for samp in samp_list:
						if samp in line:
							clear = False
							break
					
					if clear == True:
						fo.write(line)
		

###############
# matrix

class matrixApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Name for Output File Required*')


		with open(InFile, 'r') as fi:  # read in VCF file and appends all non-header lines 
			file = fi.readlines()      # to a list-- 'vcf'. Each item is a new line.

		head_bool = True 

		for line in file:
			if '##' not in line:

				if head_bool == True:
					isolates = line.split('\t')[9:] # extracts isolate names from VCF. List contains isolate names in order
					isolates[-1] = isolates[-1][:-1] # removes newline at end of string
					isolates = list(map(lambda x:[x], isolates)) 
					head = ['']
					head_bool = False

				else:
					name = str(line.split('\t')[0]) + ':' + str(line.split('\t')[1]) # extracts variants scaffold and position
					head.append(name) # unique variant name is 'scaffold':'position'

					for x in range(len(isolates)): # cycles through isolate columns list
						temp = isolates[x]
						temp.append(line.split('\t')[9+x].split(':')[0].strip())  # takes genotype assignment for isolate on each line
						isolates[x] = temp				
						

		header = '\t'.join(head) # creates tab separated string containing variant ID's.

		with open(OutFile, 'w') as fo:  # creates tab separated tables where columns are variants and rows are isolate genotypes
			fo.write(header + '\n')
			for row in isolates:
				fo.write('\t'.join(row) + '\n')


###############
# pat_match

class pat_matchApp():

	def __init__(self):
		self.verbose = False

	def pat_match(pattern, text):  #  matching algorithm, returns either '1' (match) or '0' (no match)

			pattern = pattern.split(',')  

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
			raise Exception('*VCF File or VCF Table Required*')
		elif Pattern == 'Error2':
			raise Exception('*Presence/Absence Pattern Required*')
		elif OutFile == 'Error3':
			raise Exception('*Name for Output File Required*')


		with open(InFile, 'r') as fi:
			file = fi.readlines()


		with open(OutFile, 'w') as fo:  # creates a new file with every line being a gene header + '\n'

			if '##' in file[1]: # Checks if File is in VCF or Table format

				for line in file:
					if '#' not in line:
						genotype_map = map(lambda x : x.split(':')[0].strip(), line.split('\t')[9:])
						genotype_calls = list(genotype_map)
						result = pat_matchApp.pat_match(Pattern, genotype_calls)

						if result == 1:
							name = str(line.split('\t')[0]) + ':' + str(line.split('\t')[1]) # extracts variants scaffold and position
							fo.write(name + '\n')

			else:

				l = len(file[0].split('\t'))  # length variable for loop below	


				for x in range(1, l):  # parses tab deliminated text file for pattern and generates a list of texts to search 
					temp = []
					for y in range(1, len(file)):  # splits ever line of file into a list by tabs, assigns each item on each line at a given position to a list
						number = file[y].split('\t')[x] 
						temp.append(number.strip()) 

					result = pat_matchApp.pat_match(Pattern, temp)  # outputs a '1' or '0' and appends to a list 

					if result == 1:
						fo.write(file[0].split('\t')[x].strip() + '\n')


###############
# clear

class clearApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF or Table File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Name for Output File Required*')

		var_lines =[]

		with open(InFile, 'r') as fi: # reads in VCF table file 
			file = fi.readlines() 

		with open(OutFile, 'w') as fo:  # creates file with same hash lines as input file and only variant lines without no calls

			if '##' in file[1]: # Checks if File is in VCF or Table format

				for line in file:
					if '#CHROM' in line:
						fo.write(line)  # writes all header lines to file
					else:
						for y in range(9, len(line.split('\t'))):  # iterates over every column in every variant line
							if '.' in line.split('\t')[y][0]: 
								break
						else:
							fo.write(line)

			else:

				no_calls_list = []

				header = file[0]
				file = file[1:]

				for line in file:
					line_list = line.split('\t')

					for i in range(1,len(line_list)):
						if '.' in line_list[i]:
							no_calls_list.append(i)

				file = [header] + file 
							
				with open(OutFile, 'w') as fo:  # creates file with same hash lines as input file and only variant lines without no calls
					for line in file:
						line_list = line.split('\t')
						for i in no_calls_list:
							if i != len(line_list) - 1:
								line_list = line_list[:i] + line_list[i+1:]
							else: 
								line_list = line_list[:i] 
								line_list[-1] = line_list[-1] + '\n'
						fo.write('\t'.join(line_list))


###############
# samp_comp

class samp_compApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, Samples, Exclude, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF Table File Required*')
		elif Samples == 'Error2':
			raise Exception('*Samples Required*')

		if len(Samples.split(',')) < 2:
				raise Exception('*More Than One Sample Required*')
		else:
			samples = set(Samples.split(','))

		try:	
			ex_samples = Exclude.split(',')
			print('Excluded Samples: ', end="")
			print(ex_samples)

		except:
			ex_samples = []
	

		if len(samples) == 1 and len(ex_samples) == 0:
			print('Shared Variants: ' + str(int(len(file[0].split('\t')))-1))
			print('Similarity: 1.0') 
			return

		sim = 0 
		sim_list = []


		with open(InFile, 'r') as fi: # reads in VCF table from vat_mat command 
			file = fi.readlines()


		if '##' in file[1]:  # Checks if File is in VCF or Table format

			in_indices = []
			ex_indices = []
			comp = False
			total = 0

			for line in file:

				if '#CHROM' in line:

					split_line = line.split('\t')

					comp = True

					for x in range(9, len(split_line)):

						if split_line[x].strip() in samples:
							in_indices.append(x)
						
						elif ex_samples != [] and split_line[x].strip() in ex_samples: 
							ex_indices.append(x)

					if len(in_indices) != len(samples):
						raise Exception('*One or More Comparison Samples Not Found in VCF*')

					if len(ex_indices) != len(ex_samples):
						raise Exception('*One or More Exclude Samples Not Found in VCF*')

				
				elif comp == True:

					total += 1

					split_line = line.split('\t')

					include = True

					ref_genotype = split_line[in_indices[0]].split(':')[0]

					for i in range(1, len(in_indices)):

						if ref_genotype != split_line[in_indices[i]].split(':')[0]:
							include = False
							break

					if include == True and len(ex_indices) != 0:

						for j in range(len(ex_indices)):

							if ref_genotype == split_line[ex_indices[j]].split(':')[0]:
								include = False
								break

					if include == True:

						sim += 1 
						sim_list.append(str(split_line[0]) + ":" + str(split_line[1]))

			print('Total Variants: ' + str(total))

		else:

			ex_list = []
			samp_list = []


			for i in range(len(file)):  # iterates over every row
				
				if file[i].split('\t')[0] in samples:  # finds line with corresponding first sample ID
					samp_list.append(file[i].split('\t')[1:])  # creates list with every genotype call for every variant for first isolate
					
					if len(samples) == 1:
						samp_list.append(file[i].split('\t')[1:])

				elif file[i].split('\t')[0] in ex_samples:
					ex_list.append(file[i].split('\t')[1:])
			

			if len(samp_list) != len(samples):
				raise Exception('*One or More Comparison Samples Not Found in VCF Table*')

			if len(ex_list) != len(ex_samples):
				raise Exception('*One or More Exclude Samples Not Found in VCF Table*')


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

			print('Total Variants: ' + str(int(len(file[0].split('\t')))-1))  # counts number of columns 

		if Exclude == None:
			print('Shared Variants: ' + str(sim)) 
		else:
			print('Variants not shared with excluded samples: ' + str(sim))

		if OutFile != None:
			with open(OutFile, "w") as fo:
				for sim in sim_list:
					fo.write(sim + '\n')


###############
# sim_mat

class sim_matApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF or Table File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Output File Name Required*')

		table = {}
		samples = []
		sim_table = []

		with open(InFile, "r") as fi:  # reads in vcf matrix
			file = fi.readlines()


		if '##' in file[1]: # Checks if File is in VCF or Table format

			for i in range(len(file)):  # creates dictionary with sample as key an allele sequence as value 	

				if "#CHROM" in file[i]:
					pos = i
					for sample in file[i].split("\t")[9:]:
						table[sample.strip()] = [0] * len(file[i].split("\t")[9:])
						samples.append(sample.strip())

					break


			num_samps = len(samples)
			var_length = len(file[pos+1:])

			for i in range(num_samps):
				sample = samples[i]
				for j in range(pos+1, len(file)):
					for k in range(num_samps):
						if file[j].split("\t")[9+i].split(':')[0] == file[j].split("\t")[9+k].split(':')[0]:
							table[sample][k] += 1

			
			for samp in samples:  # populates the distance matrix with 0's to the propper dimensions
				sim_table.append(list(map(lambda x : x / var_length, table[samp])))

			with open(OutFile, "w") as fo:  # opens output file
				fo.write(' ' + '\t' + '\t'.join(samples) + '\n')  # populates matrix column headers
				for x in range(len(sim_table)):           # populates matrix by row
					temp = samples[x] 
					for y in range(len(sim_table[x])):
						temp += '\t' + str(sim_table[x][y])
					fo.write(temp + '\n') 


		else:

			for line in file:  # creates dictionary with sample as key an allele sequence as value 
				table[line.split("\t")[0]] = line.split("\t")[1:]
				samples.append(line.split("\t")[0])  # also adds all samples to a list

			num_samps = len(samples) 

			for i in range(num_samps - 1):  # populates the distance matrix with 0's to the propper dimensions
				sim_table.append([])
				for j in range(num_samps - 1):
					sim_table[i].append("0")

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
					sim_table[x-1][y-1] = perc_sim  # adds above value to correct spot in matrix

			with open(OutFile, "w") as fo:  # opens output file
				fo.write('\t'.join(samples) + '\n')  # populates matrix column headers
				for x in range(len(sim_table)):           # populates matrix by row
					temp = samples[x+1] 
					for y in range(len(sim_table[x])):
						temp += '\t' + str(sim_table[x][y])
					fo.write(temp + '\n') 

###############
# allele_seq

class allele_seqApp():

	def __init__(self):
		self.verbose = False


	# removes pesky quotes from wieirdly formatted VCFs
	def unquote(string): 

		if "\"" in string:
			string = string.replace('\"', '')

		if "\'" in string:
			string = string.replace('\'', '')

		return string 

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Output File Name Required*')

		with open(InFile, "r") as fi:  # reads in vcf matrix
			file = fi.readlines()


			line_num = 0

			for line in file:

				if "#CHROM" in line:

					samples = line.split("\t")[9:] # samples list
					sequences = [""] * len(samples) # empty sequence list

				elif "#" not in line and len(samples) != 0:

					temp = line.split("\t") # columns 

					for i in range(9, len(temp)): # iterates over sample columns


						geno_field = allele_seqApp.unquote(temp[i])
						allele = geno_field.split(":")[0] 

						if allele == "0":

							sequences[i - 9] += allele_seqApp.unquote(temp[3])

						elif allele == "." or "." in allele:

							sequences[i - 9] += "-"

						else:

							alt = int(allele) - 1
							sequences[i - 9] += allele_seqApp.unquote(temp[4].split(',')[alt])


			with open(OutFile, "w") as fo: # writes to output multifasta

				for i in range(len(samples)):
							
					fo.write(">" + allele_seqApp.unquote(samples[i].strip()) + "\n")
					fo.write(sequences[i].strip() + "\n")


#####################################################################################

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
   		app = removeApp()
   		return app.start(args.InFile, args.Sample, args.OutFile)

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
   		app = matrixApp()
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
   		app = pat_matchApp()
   		return app.start(args.InFile, args.Pattern, args.OutFile)

###############
# clear
 
def clearParser(subparsers):
	clear_parser = subparsers.add_parser('clear',
		help='Removes all variants from VCF files or a Variant Matrix where a genotype call could not be made for any isolate')
	clear_parser.add_argument('-i', '--input', help='VCF or Table file', dest='InFile', type=str, default='Error1')
	clear_parser.add_argument('-o', '--output', help='Name of output table', dest='OutFile', type=str, default='Error2')
	
	return clear_parser

class clearCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = clearApp()
   		return app.start(args.InFile, args.OutFile)

###############
# samp_comp
  
def samp_compParser(subparsers):
	samp_comp_parser = subparsers.add_parser('samp_comp',
		help='Finds the percent of variants shared between two samples in a VCF file or Variant Matrix')
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
# sim_mat
  
def sim_matParser(subparsers):
	sim_mat_parser = subparsers.add_parser('sim_mat',
		help='Create a similarity matrix of all samples in a a VCF file or Variant Matrix')
	sim_mat_parser.add_argument('-i', '--input', help='VCF or Table File(output from vcf_mat)', dest='InFile', type=str, default='Error1')
	sim_mat_parser.add_argument('-o', '--output', help='Name of output similarity matrix', dest='OutFile', type=str, default='Error2')

	return sim_mat_parser

class sim_matCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = sim_matApp()
   		return app.start(args.InFile, args.OutFile)

###############
# allele_seq
  
def allele_seqParser(subparsers):
	allele_seq_parser = subparsers.add_parser('allele_seq',
		help='Creates a multi-fasta file containing the sequence of all snp alleles for each sample in a VCF (haploid and snps only for now)')
	allele_seq_parser.add_argument('-i', '--input', help='VCF file', dest='InFile', type=str, default='Error1')
	allele_seq_parser.add_argument('-o', '--output', help='Name of output file', dest='OutFile', type=str, default='Error2')

	return allele_seq_parser

class allele_seqCMD():

	def __init__(self):
		pass

	def execute(self, args):
   		app = allele_seqApp()
   		return app.start(args.InFile, args.OutFile)

#####################################################################################

###############
# argument parser

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Miscellaneous functions for manipulating and analyzing VCF files.', add_help=True,
        epilog="For questions or comments, contact Bradley Jenner <bnjenner@ucdavis.edu>")
    subparsers = parser.add_subparsers(help='commands', dest='command')
    removeParser(subparsers)
    matrixParser(subparsers)
    pat_matchParser(subparsers)
    clearParser(subparsers)
    samp_compParser(subparsers)
    sim_matParser(subparsers)
    allele_seqParser(subparsers)
    args = parser.parse_args()
    return args

def main():
	remove = removeCMD()
	matrix = matrixCMD()
	pat_match = pat_matchCMD()
	clear = clearCMD()
	samp_comp = samp_compCMD()
	sim_mat = sim_matCMD()
	allele_seq = allele_seqCMD()
	commands = {'remove': remove, 'matrix': matrix, # remove, matrix
				'pat-match': pat_match, 'clear': clear, # pat-match, clear
				'samp_comp': samp_comp , 'sim_mat': sim_mat, # compare, sim-matrix
				'allele_seq': allele_seq} # allele-seq
	args = parseArgs()
	commands[args.command].execute(args)
	
if __name__ == '__main__':
	main()
