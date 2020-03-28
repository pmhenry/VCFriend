from . import utilities

###############
# compare

class compareApp():

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


		with open(InFile, 'r') as fi: # reads in VCF / matrix
			file = fi.readlines()


		if '##' in file[1]:  # Checks if File is in VCF or Table format

			in_indices = []
			ex_indices = []
			comp = False
			total = 0

			for line in file:

				if '#CHROM' in line:

					split_line = utilities.unquote(line).split('\t')

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

					split_line = utilities.unquote(line).split('\t')

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

				line = utilities.unquote(file[i])
				split_line = line.split('\t')[0] 

				if split_line[0] in samples:  # finds line with corresponding first sample ID
					samp_list.append(split_line[1:])  # creates list with every genotype call for every variant for first isolate
					
					if len(samples) == 1:
						samp_list.append(line.split('\t')[1:])

				elif split_line[0] in ex_samples:
					ex_list.append(split_line[1:])
			

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
