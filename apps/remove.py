from . import utilities

###############
# remove

class removeApp():

	def __init__(self):
		self.verbose = False


	def start(self, InFile, Samples, OutFile):
	
		if InFile == 'Error1':
			raise Exception('*VCF or Table File Required*')
		elif Samples == 'Error2':
			raise Exception('*Sample IDs Required*')
		elif OutFile == 'Error3':
			raise Exception('*Name for Output File Required*')

		var_lines = []
		samp_list = []


		if ".txt" in Samples:

			with open(Samples, 'r') as fi:
				samp_list = fi.readlines()
				samp_list = list(map(lambda x : x.strip(), samp_list))
		else:

			samp_list = Samples.split(',')


		if len(samp_list) == 0:
			raise Exception('*Sample IDs Required*')

		with open(InFile, 'r') as fi: # reads in VCF / matrix
			file = fi.readlines()


		with open(OutFile, 'w') as fo:  # creates file with original hash lines from input VCF.

			if '##' in file[1]: # Checks if File is in VCF or Table format

				for line in file:

					if '##' in line:
						fo.write(line)
					
					else:
					
						if '#CHROM' in line:

							header = utilities.unquote(line).split('\t')
							header = list(map(lambda x : x.strip(), header))
							rm_list = []
					
							for x in range(len(header) - 1, 8, -1):

								samp = header[x].strip()
								if samp in samp_list:
									rm_1 = x # assigns index to remove in future lines
									rm_list.append(rm_1)
									header.remove(header[x]) # removes column corresponding to removed sample

							fo.write('\t'.join(header) + "\n") 
					
						else:

							temp = utilities.unquote(line).split('\t')  # creates temp list where each item is a column from the VCF file
							temp  = list(map(lambda x : x.strip(), temp))

							for col in range(len(temp) -1, -1, -1):

								if col in rm_list:

									temp.remove(temp[col])  # removes column corresponding to removed sample
							
							fo.write('\t'.join(temp) + "\n")

			else:

				for line in file:
					clear = True
					
					for samp in samp_list:
					
						if samp in line:
							clear = False
							break
					
					if clear == True:
						fo.write(utilities.unquote(line))
		
