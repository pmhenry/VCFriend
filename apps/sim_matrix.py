from . import utilities

###############
# sim-matrix

class sim_matrixApp():

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

		with open(InFile, "r") as fi:  # reads in vcf / matrix
			file = fi.readlines()


		if '##' in file[1]: # Checks if File is in VCF or Table format

			for i in range(len(file)):  # creates dictionary with sample as key an allele sequence as value 

				if "#CHROM" in file[i]:
					line = utilities.unquote(file[i])
					pos = i

					for sample in line.split("\t")[9:]:
						table[sample.strip()] = [0] * len(line.split("\t")[9:])
						samples.append(sample.strip())

					break


			num_samps = len(samples)
			var_length = len(file[pos+1:])

			for i in range(num_samps):
				sample = samples[i]

				for j in range(pos+1, len(file)):
					
					for k in range(num_samps):

						if utilities.unquote(file[j]).split("\t")[9+i].split(':')[0] == utilities.unquote(file[j]).split("\t")[9+k].split(':')[0]:
							table[sample][k] += 1

			
			for samp in samples:  # populates the distance matrix with 0's to the propper dimensions
				sim_table.append(list(map(lambda x : x / var_length, table[samp])))

			with open(OutFile, "w") as fo:  # opens output file
				fo.write(' ' + '\t' + utilities.unquote('\t'.join(samples)) + '\n')  # populates matrix column headers
				
				for x in range(len(sim_table)):           # populates matrix by row
					temp = samples[x] 
					
					for y in range(len(sim_table[x])):
						temp += '\t' + str(sim_table[x][y])
					
					fo.write(temp + '\n') 


		else:

			for line in file:  # creates dictionary with sample as key an allele sequence as value 

				temp = utilities.unquote(line)
				temp_split = temp.split("\t")


				table[temp_split[0]] = temp_split[1:]
				samples.append(temp_split[0])  # also adds all samples to a list

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

