import sys
import os
import re


project_name = config['project_name']
workflow_path = config['workflow_path']
input_files = config['input_files']
chromosomes = config['chromosomes']
accession_mapping_file = config['accession_mapping_file']
accession_key = config['accession_key']
phenotype_key = config['phenotype_key']
phenotype_file = config['phenotype_file']
phenotypes = config['phenotypes']
p_value_filtering_threshold = config['p_value_filtering_threshold']
gff_file = config['gff_file']
gff_category = config['gff_category']
gff_key = config['gff_key']
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
    possible_sample = re.sub('(\\.vcf.*)','',str(os.path.basename(input_files[i])))
    if possible_sample not in input_samples:
        input_samples.append(possible_sample)
    possible_extension = re.sub(possible_sample,'',str(os.path.basename(input_files[i])))
    if possible_extension != input_extension:
        input_extension = possible_extension
    possible_sample = re.sub(re.compile('_' + str(chromosomes[i]) + '$'),'',possible_sample)
    if input_sample != possible_sample:
        input_sample = possible_sample


rule all:
    input:
        expand(os.path.join(os.path.abspath(output_folder),'grep_effects',input_sample + '_{chromosome}' + input_extension),chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder),'generate_phenotype_distribution_with_statistical_testing',input_sample + '_{chromosome}__{phenotype}.txt'),chromosome=chromosomes,phenotype=phenotypes),
        expand(os.path.join(os.path.abspath(output_folder),'combine_phenotype_distribution_statistical_results_by_chromosomes','{project_name}_{{chromosome}}.txt'.format(project_name=project_name)),chromosome=chromosomes),
        expand(os.path.join(os.path.abspath(output_folder),'combine_phenotype_distribution_statistical_results_by_phenotypes','{project_name}_{{phenotype}}.txt'.format(project_name=project_name)),phenotype=phenotypes),
        os.path.join(os.path.abspath(output_folder),'summarize_phenotype_distribution_statistical_results','{project_name}.txt'.format(project_name=project_name))


rule grep_input_effects:
    input:
        in_file=os.path.join(os.path.abspath(input_folder),input_sample + '_{chromosome}' + input_extension)
    output:
        out_file=os.path.join(os.path.abspath(output_folder),'grep_effects',input_sample + '_{chromosome}' + input_extension)
    shell:
        """
        zgrep -e "^#" -e "^#CHROM" -e "frameshift_variant" -e "exon_loss_variant" -e "duplication" -e "inversion" \
        -e "feature_ablation" -e "gene_fusion" -e "rearranged_at_DNA_level" -e "missense_variant" -e "protein_protein_contact" \
        -e "structural_interaction_variant" -e "rare_amino_acid_variant" -e "splice_acceptor_variant" -e "splice_donor_variant" \
        -e "stop_lost" -e "start_lost" -e "stop_gained" -e "inframe_insertion" -e "disruptive_inframe_insertion" -e "inframe_deletion" \
        -e "disruptive_inframe_deletion" {input.in_file} | grep -v -e "DEL" -e "INS" | bgzip > {output.out_file}
        """


include: '../rules/python/generate_phenotype_distribution_with_statistical_testing.smk'
include: '../rules/python/combine_phenotype_distribution_statistical_results_by_chromosomes.smk'
include: '../rules/python/combine_phenotype_distribution_statistical_results_by_phenotypes.smk'
include: '../rules/python/summarize_phenotype_distribution_statistical_results.smk'
