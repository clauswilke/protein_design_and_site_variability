Written/Last Updated on October 7, 2013 by Eleisha Jackson
In this directory contains all of the scripts needed for the "Amino-acid site variability among natural and designed proteins" Paper
 
Directories:
alignment_scripts/
	Provides the scripts that were used align the yeast-proteins data set 
duncan_analysis/
	Contains all of the scripts used to analyze the yeast-proteins data set
	analysis_functions.py
		This a helper file that contains a series of functions used in the analysis				
	calculate_distribution_files.py
		This is a script that calculates all of the amino acid count data for the designed yeast proteins.
	calculate_natural_distribution_files.py
		This is a script that calculates all of the amino acid count data for the natural yeast proteins.
	generate_graph_data.py
		This is a script that calculates all of the entropy, KL-Divergence, and RSA-Entropy correlation analysis
	unordered_analysis/
		Contains the scripts that are used to analyze the yeast-proteins data set (This is for the unranked data)
	ordered_analysis/
		Contains the scripts that are used to analyze the yeast-proteins data set (This is for the ranked data)

graphing/
	A series of python scripts that analysis the data and produce the figures for the paper.
	analysis_functions.py
		This a helper file that contains a series of functions used in the analysis				
	combo_mixed_analysis.py
		This script is the script that is used to create Figure 6.
	get_combo_freq_plots_noah.py
		This script is the script that is used to create Figures S3, S4, and S5.
	get_combo_freq_plots.py
		This script is the script that is used to create Figures 2, S1, and S2.			
	get_mean_data_boxplots_noah.py
		This script is the script that is used to create Figure S6, and S8.
	get_mean_data_boxplots.py
		This script is the script that is used to create Figure 3, and 4.
	get_mean_entropy_combo_plot.py
		This script is the script that is used to create Figure 1.
	get_mean_entropy_median_plot_combo.py
		This script is the script that is used to create Figure 5.
	get_RSA_entropy_Combination_Plot.py
		This script is the script that is used to create Figure S7.

noah_analysis/
	Contains all of the scripts used to analyze the protein-domains dataset
	analysis_functions.py
		This a helper file that contains a series of functions used in the analysis.				
	calculate_distribution_files_noah.py
		This is a script that calculates all of the amino acid count data for the designed protein domains.
	calculate_distribution_files_soft_noah.py
		This is a script that calculates all of the amino acid count data for the "soft" designed protein domains.
	calculate_natural_distribution_files_noah.py
		This is a script that calculates all of the amino acid count data for the natural protein domains.
	generate_graph_data.py
		This is a script that calculates all of the entropy, KL-Divergence, and RSA-Entropy correlation analysis.
	ordered_analysis/
		These are the scripts used in the ranked analysis 
	unordered_analysis/
		These are the scripts used in the unranked analysis

r_scripts/
	Contains all of the R scripts used in the analysis
	calculate_stats.R 
		This is a script that performs the statistics in the manuscript. 
	The following files are needed to calculate the statistics:
	Files that are created from the generate_graph_data.py file: 
		mixed_duncan_stats_file.csv
		stats_data.csv 
	Files that are created from the generate_graph_data_noah.py file:
		mixed_noah_stats_file.csv
		stats_data_noah.csv

rosetta_scripts/
	Contains the two shell scripts that were used to run RosettaDesign and Backrub
	ucsf_yeash.sh
		This is a shell script that was used to create the flexible backbon protein designs
		For arguments you must give it the PDB name and the temperature used for Backrub
	ucsf_yeash_fixed.sh
		This the shell script that was used to create the fixed backbone protein design
		For arguments you must give it the PDB name and the temperature used for Backrub.
		Since it is a fixed design we just set the temperature to 0.0 purely for identification purposes.

sequences/
	Contains all of the sequences (both designed and natural) for the project
	duncan_sequences/	
		aligned_sequences/
			These are the natural alignments for OUR dataset used in the analysis
		designed_sequences/
			These are the designed sequences extracted from the yeast-protein designed proteins 
	noah_sequences/
		backrub_sequences/
			These are the designed sequences extracted from the protein-domain designed proteins
		natural_alignments/
			These are the natural alignments for OUR dataset used in the analysis

structures/
	Contains all of the protein structures used for the design process
	duncan_structures/
		This contains the PDB Files that were used as templates in the yeast-proteins design process
	noah_structures/
		This contains the PDB Files that were used as templates in the protein-domain design process


	