rule summarize_phenotype_distribution_statistical_results:
    input:
        in_file=expand(os.path.join(os.path.abspath(output_folder),'combine_phenotype_distribution_statistical_results_by_chromosomes','{project_name}_{{chromosome}}.txt'.format(project_name=project_name)),chromosome=chromosomes)
    params:
        all_in_files=' '.join(['-i ' + os.path.join(os.path.abspath(output_folder),'combine_phenotype_distribution_statistical_results_by_chromosomes','{project_name}_{chromosome}.txt'.format(project_name=project_name,chromosome=chromosome)) for chromosome in chromosomes])
    output:
        out_file=os.path.join(os.path.abspath(output_folder),'summarize_phenotype_distribution_statistical_results','{project_name}.txt'.format(project_name=project_name))
    log:
        os.path.join(os.path.abspath(output_folder),'summarize_phenotype_distribution_statistical_results_log','{project_name}.log'.format(project_name=project_name))
    threads: threads
    resources:
        memory=memory
    shell:
        """
        python3 {workflow_path}/scripts/python/summarize_phenotype_distribution_statistical_testing_results.py \
        {params.all_in_files} \
        -o {output.out_file} 2>&1 > {log}
        """
