from . import utilities

###############
# pat_match

class pat_matchApp():

	def __init__(self):
		self.verbose = False


	def start(self, InFile, Pattern, OutFile):


		if InFile == 'Error1':
			raise Exception('*VCF File or VCF Table Required*')
		elif Pattern == 'Error2':
			raise Exception('*Presence/Absence Pattern Required*')
		elif OutFile == 'Error3':
			raise Exception('*Name for Output File Required*')


		with open(InFile, 'r') as fi: # reads in VCF / matrix
			file = fi.readlines()


		with open(OutFile, 'w') as fo:  # creates a new file with every line being a gene header + '\n'

			if '##' in file[1]: # Checks if File is in VCF or Table format

				for line in file:

					if '#' not in line:
						genotype_map = map(lambda x : x.split(':')[0].strip(), utilities.unquote(line).split('\t')[9:])
						genotype_calls = list(genotype_map)
						result = utilities.pat_match(Pattern, genotype_calls)

						if result == 1:
							name = utilities.unquote(str(line.split('\t')[0]) + ':' + str(line.split('\t')[1])) # extracts variants scaffold and position
							fo.write(name + '\n')

			else:

				l = len(file[0].split('\t'))  # length variable for loop below	


				for x in range(1, l):  # parses tab deliminated text file for pattern and generates a list of texts to search 
					temp = []
					
					for y in range(1, len(file)):  # splits ever line of file into a list by tabs, assigns each item on each line at a given position to a list
						number = utilities.unquote(file[y]).split('\t')[x] 
						temp.append(number.strip()) 

					result = utilities.pat_match(Pattern, temp)  # outputs a '1' or '0' and appends to a list 

					if result == 1:
						fo.write(file[0].split('\t')[x].strip() + '\n')

