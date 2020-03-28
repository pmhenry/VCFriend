from . import utilities

###############
# allele_seq

class allele_seqApp():

	def __init__(self):
		self.verbose = False


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

					samples = utilities.unquote(line).split("\t")[9:] # samples list
					sequences = [""] * len(samples) # empty sequence list

				elif "#" not in line and len(samples) != 0:

					temp = line.split("\t") # columns 

					for i in range(9, len(temp)): # iterates over sample columns


						geno_field = utilities.unquote(temp[i])
						allele = geno_field.split(":")[0] 

						if allele == "0":

							sequences[i - 9] += utilities.unquote(temp[3])

						elif allele == "." or "." in allele:

							sequences[i - 9] += "-"

						else:

							alt = int(allele) - 1
							sequences[i - 9] += utilities.unquote(temp[4]).split(',')[alt]


			with open(OutFile, "w") as fo: # writes to output multifasta

				for i in range(len(samples)):
							
					fo.write(">" + samples[i].strip() + "\n")
					fo.write(sequences[i].strip() + "\n")

