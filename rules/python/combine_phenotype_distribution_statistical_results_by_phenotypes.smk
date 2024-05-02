rule combine_phenotype_distribution_statistical_results_by_phenotypes:
    input:
        in_file=expand(os.path.join(os.path.abspath(output_folder),'generate_phenotype_distribution_with_statistical_testing',input_sample + '_{chromosome}__{{phenotype}}.txt'),chromosome=chromosomes)
    params:
        all_in_files=' '.join(['-i ' + os.path.join(os.path.abspath(output_folder),'generate_phenotype_distribution_with_statistical_testing',input_sample + '_{chromosome}__{{phenotype}}.txt'.format(chromosome=chromosome)) for chromosome in chromosomes]),
        phenotype='{phenotype}'
    output:
        out_file=os.path.join(os.path.abspath(output_folder),'combine_phenotype_distribution_statistical_results_by_phenotypes','{project_name}_{{phenotype}}.txt'.format(project_name=project_name))
    log:
        os.path.join(os.path.abspath(output_folder),'combine_phenotype_distribution_statistical_results_by_phenotypes_log','{project_name}_{{phenotype}}.log'.format(project_name=project_name))
    threads: threads
    resources:
        memory=memory
    shell:
        """
        python3 {workflow_path}/scripts/python/combine_phenotype_distribution_statistical_testing_results.py \
        {params.all_in_files} \
        -o {output.out_file} 2>&1 > {log}
        """
