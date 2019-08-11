# Halteria_genome_research
##Bioinformatics tools used in Halteria genome research.  
##Usage of all python (or perl) scripts:  
###(1)Code S1. simulation_binomial_sampling.py  
This is a python script.  
Required module: No additional python module is required.  
Usage: python simulation_binomial_sampling.py \<input file: table recording pool population data\> \<output file\> \<total reads number\>   
  Eg. python simulation_biomial_sampling.py pol_data.tab output_data.tab 186978547  
  
The input file is a file in table format. Format is shown below.  

chromosome	position	minor_allele	MAF	major_allele	MAAF	position_type	poly_type	reads_num	gene
NODE_5302_length_4220_cov_357.479	30	T	0.264044943820225	A	0.69873595505618	inter_genic	bi	1424	NULL
NODE_5302_length_4220_cov_357.479	31	G	0.262829403606103	A	0.701109570041609	inter_genic	bi	1442	NULL
NODE_5302_length_4220_cov_357.479	34	T	0.277894736842105	A	0.72140350877193	inter_genic	bi	1425	NULL
NODE_5302_length_4220_cov_357.479	37	A	0.291666666666667	G	0.690277777777778	inter_genic	bi	1440	NULL
