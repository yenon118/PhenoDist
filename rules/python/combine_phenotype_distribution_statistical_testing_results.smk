rule combine_phenotype_distribution_statistical_testing_results:
    input:
        in_file = expand(os.path.join(os.path.abspath(output_folder), 'generate_phenotype_distribution_with_statistical_testing', input_sample+'_{{chromosome}}__{phenotype}.txt'), phenotype=phenotypes)
    params:
        all_in_files = ' '.join(['-i '+os.path.join(os.path.abspath(output_folder), 'generate_phenotype_distribution_with_statistical_testing', input_sample+'_{{chromosome}}__{phenotype}.txt'.format(phenotype=phenotype)) for phenotype in phenotypes]),
        chromosome = '{chromosome}'
    output:
        out_file = os.path.join(os.path.abspath(output_folder), 'combine_phenotype_distribution_statistical_testing_results', '{project_name}_{{chromosome}}.txt'.format(project_name=project_name))
    log:
        os.path.join(os.path.abspath(output_folder), 'combine_phenotype_distribution_statistical_testing_results_log', '{project_name}_{{chromosome}}.log'.format(project_name=project_name))
    threads: threads
    resources:
        memory = memory
    shell:
        """
        python3 {workflow_path}/scripts/python/combine_phenotype_distribution_statistical_testing_results.py \
        {params.all_in_files} \
        -o {output.out_file} 2>&1 > {log}
        """