###############
# clean

class cleanApp():

	def __init__(self):
		self.verbose = False

	def start(self, InFile, OutFile):

		if InFile == 'Error1':
			raise Exception('*VCF or Table File Required*')
		elif OutFile == 'Error2':
			raise Exception('*Name for Output File Required*')

		var_lines =[]

		with open(InFile, 'r') as fi: # reads in VCF / matrix
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
