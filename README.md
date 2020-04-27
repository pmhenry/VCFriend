# VCFriend

Variant calling softwares utilizes VCF files that can sometimes be difficult to parse. VCFriend is a set of tools that can make extracting and analyzing information from VCF files easier.
	
## Installation:

To install, run following commands:

	git clone https://github.com/bnjenner/VCFriend.git
	echo "path/to/repo/" >> ~/.bash_profile # or equivalent 
	source ~/.bash_profile
	VCFriend.py -h # help page means installation was successful.
---
## Usage 

VCFriend.py [FUNCTION] [OPTIONS]

## Functions

  allele-seq:

	Creates a multi-fasta file containing the sequence of all snp alleles for each sample in a VCF (haploid and snps only for now).

	OPTIONS:
	-i [ --input ]          Input VCF
	
	-o [ --output ]         Name of Output Fasta File


  clean:
  
	Removes variants for which a genotype could not be made for all samples in VCF or Table file. Returns new VCF/Matrix file without specified sample.

	OPTIONS:
	-i [ --input ] 		Input VCF or Matrix

	-o [ --output ] 	Name of Output File		


  compare:
  
	Returns the number of variants shared between any number of samples excluding 0 or more samples in VCF files or variant matrix (vcf_mat output). Function is phase sensitive for genotypes.

	OPTIONS:
	-i [ --input ] 		Input VCF File or Matrix

	-s [ --samples ]	List of samples to compare separated by commas 

	-x [ --exclude ]	List of samples to exclude from comparison separated by commas

	-o [ --output ]         Output list of shared variant names


  matrix:
  
	VCF file simplification. Generates variant matrix (tab delimited text file) containing only genotype data. Samples are rows and variants are columns. 

	VCF Table Example: 

		  	  variation_1   variation_2     variation_3    variation_4
    
    	genome_1            1             0               1              1
 
   		genome_2            1             1               0              1

    	genome_3            0             1               1              1

    	genome_4            0             0               0              0
 
    	genome_5            1             0               1              0
 
    	genome_6            1             0               1              1


	OPTIONS:
	-i [ --input ] 		Input VCF File

	-o [ --output ] 	Name of Output File


  pat-match:
  
	Searches VCF file or variant matrix (matrix output) for patterns of presence/absence for variants among genomes. Returns file with list of variants matching the designated pattern ofpresence/absence. For patterns, 1's indicate that variant must be present in given genoms while 0's indicate absence. N's exclude the genome from the pattern match. N's in pattern mean sample is ignored. Y's require that a genotype call is made for the sample. .'s correspond to no genotype calls in the VCF.

  	OPTIONS:
	-i [ --input ] 		Input VCF File or Matrix

	-p [ --pattern ]	Comma separated string in quotes containing genotype calls or N's, Y's, and .'s. Pattern is phase sensitie. Sequence of search pattern must be correspond to order of samples in VCF Table

	-o [ --output ] 	Name of Output File


   remove:
   
	Removes samples from a VCF or Variant Matrix file.	
  	
	OPTIONS:
	-i [ --input ] 		Input VCF or Matrix File

	-s [ --sample ]		Comma separated list of ID of Samples to be removed or text file with a sample on each line

	-o [ --output ] 	Name of Output File


  sim-matrix:
  
	Returns similarity matrix for samples.

	OPTIONS:
	-i [ --input ]          Input VCF or Matrix

	-o [ --output ]         Name of Output File
  
---

Last Updated: March 26th, 2020

For comments of questions, please contact Bradley Jenner <bnjenner@ucdavis.edu>
