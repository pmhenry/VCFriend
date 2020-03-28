from . import utilities

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
					isolates = utilities.unquote(line).split('\t')[9:] # extracts isolate names from VCF. List contains isolate names in order
					isolates[-1] = isolates[-1][:-1] # removes newline at end of string
					isolates = list(map(lambda x:[x], isolates)) 
					head = ['']
					head_bool = False

				else:
					name = utilities.unquote(str(line.split('\t')[0]) + ':' + str(line.split('\t')[1])) # extracts variants scaffold and position
					head.append(name) # unique variant name is 'scaffold':'position'

					for x in range(len(isolates)): # cycles through isolate columns list
						temp = isolates[x]
						temp.append(utilities.unquote(line).split('\t')[9+x].split(':')[0].strip())  # takes genotype assignment for isolate on each line
						isolates[x] = temp				
						

		header = utilities.unquote('\t'.join(head)) # creates tab separated string containing variant ID's.

		with open(OutFile, 'w') as fo:  # creates tab separated tables where columns are variants and rows are isolate genotypes
			fo.write(header + '\n')
			for row in isolates:
				fo.write('\t'.join(row) + '\n')
				
