from . import utilities

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
					
							for x in range(9, len(header)):
					
								if header[x] in samp_list:
									rm_1 = x # assigns index to remove in future lines
									header.remove(header[rm_1]) # removes column corresponding to removed sample
									fo.write('\t'.join(header))
									break
					
						else:
							temp = utilities.unquote(line).split('\t')  # creates temp list where each item is a column from the VCF file
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
						fo.write(utilities.unquote(line))
		
