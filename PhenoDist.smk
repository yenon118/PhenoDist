import sys
import os
import re


project_name = config['project_name']
workflow_path = config['workflow_path']
input_files = config['input_files']
chromosomes = config['chromosomes']
beagle_window = config['beagle_window']
genome_version = config['genome_version']
accession_mapping_file = config['accession_mapping_file']
accession_key = config['accession_key']
phenotype_key = config['phenotype_key']
phenotype_file = config['phenotype_file']
phenotypes = config['phenotypes']
p_value_filtering_threshold = config['p_value_filtering_threshold']
output_folder = config['output_folder']
memory = config['memory']
threads = config['threads']


input_folder = ''
input_sample = ''
input_samples = []
input_extension = ''


for i in range(len(input_files)):
	if os.path.dirname(input_files[i]) != input_folder:
		input_folder = os.path.dirname(input_files[i])
	possible_sample = re.sub('(\\.vcf.*)', '', str(os.path.basename(input_files[i])))
	if possible_sample not in input_samples:
		input_samples.append(possible_sample)
	possible_extension = re.sub(possible_sample,'',str(os.path.basename(input_files[i])))
	if possible_extension != input_extension:
		input_extension = possible_extension
	possible_sample = re.sub(re.compile('_'+str(chromosomes[i])+'$'), '', possible_sample)
	if input_sample != possible_sample:
		input_sample = possible_sample


rule all:
	input:
		expand(os.path.join(os.path.abspath(output_folder), 'beagle_impute_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}.html'), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'snpEff_input_file', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'grep_effects', input_sample+'_{chromosome}'+input_extension), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'generate_phenotype_distribution_with_statistical_testing', input_sample+'_{chromosome}__{phenotype}.txt'), chromosome=chromosomes, phenotype=phenotypes),
		expand(os.path.join(os.path.abspath(output_folder), 'combine_phenotype_distribution_statistical_results_by_chromosomes', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name)), chromosome=chromosomes),
		expand(os.path.join(os.path.abspath(output_folder), 'combine_phenotype_distribution_statistical_results_by_phenotypes', '{project_name}_{{phenotype}}.txt'.format(project_name=project_name)), phenotype=phenotypes)


include: './rules/java/beagle_impute_input_file.smk'
include: './rules/java/snpEff_input_file.smk'

include: './rules/bash/grep_input_effects.smk'

include: './rules/python/generate_phenotype_distribution_with_statistical_testing.smk'
include: './rules/python/combine_phenotype_distribution_statistical_results_by_chromosomes.smk'
include: './rules/python/combine_phenotype_distribution_statistical_results_by_phenotypes.smk'
