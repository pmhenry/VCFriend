####################
VCFriend.py 
####################

DESCRIPTION:
	Variant calling softwares utilize VCF files that can sometimes be difficult to parse. VCFriend is a set of tools that can make analyzing and extracting information from VCF files easier.
	
USAGE:
	python3 VCFriend.py [FUNCTION] [OPTIONS]
	
FUNCTIONS: remove, vcf_mat, pat_match, nono_calls, samp_comp, dist_mat, and snp_stat
   
   REMOVE:
	Removes a sample from a VCF file.	
  	
	OPTIONS:
	-i [ --input ] 		Input VCF File

	-s [ --sample ]		ID of Sample to be removed

	-o [ --output ] 	Name of Output File
	
  VCF_MAT:
	Generates variant matrix (tab delimited text file) where samples are rows and variants are columns. Number represents genotype designation from VCF file. 

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

  PAT_MATCH:
	Searches variant matrix (vcf_mat output) for patterns of presence/absence for variants among genomes. Returns file with list of variants matching the designated pattern ofpresence/absence. For patterns, 1's indicate that variant must be present in given genoms while 0's indicate absence. N's exclude the genome from the pattern match. N's in pattern mean sample is ignored. Y's require that a genotype call is made for the sample. .'s correspond to no genotype calls in the VCF.

  	OPTIONS:
	-i [ --input ] 		Input VCF File or Table

	-p [ --pattern ]	Comma separated string in quotes containing genotype calls or N's, Y's, and .'s. Pattern is phase sensitie. Sequence of search pattern must be correspond to order of samples in VCF Table

	-o [ --output ] 	Name of Output File

  NONO_CALLS:
	Removes variants for which a genotype could not be made for all samples in VCF or Table file. Returns new VCF/Table file without specified sample.

	OPTIONS:
	-i [ --input ] 		Input VCF or Table

	-o [ --output ] 	Name of Output File		

  SAMP_COMP:
	Returns the number of variants shared between any number of samples in variant matrix (vcf_mat output).

	OPTIONS:
	-i [ --input ] 		Input VCF Table

	-s [ --samples ]	List of samples to compare separated by commas 

	-x [ --exclude ]	List of samples to exclude from comparison separated by commas

	-o [ --output ]         Output list of shared variant names

  DIST_MAT:
	Returns distance matrix from samples in variant matrix (vcf_mat output).

	OPTIONS:
	 -i [ --input ]          Input VCF Table

	 -o [ --output ]         Name of Output File
  
  SNP_STAT:
	Returns two summary table of stats at the individual variant level and averages of the same stats at the chromosome level.

	OPTIONS:
	 -i [ --input ]          Input VCF File

	 -o [ --output ]         Name of Output File Suffix
  

Last Updated: November 7th, 2019

For comments of questions, please contact Bradley Jenner <bnjenner@ucdavis.edu>
