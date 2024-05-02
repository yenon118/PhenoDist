rule generate_phenotype_distribution_with_statistical_testing:
    input:
        in_file=os.path.join(os.path.abspath(output_folder),'grep_effects',input_sample + '_{chromosome}' + input_extension),
        accession_mapping_file=accession_mapping_file,
        phenotype_file=phenotype_file
    params:
        accession_key=accession_key,
        phenotype_key=phenotype_key,
        phenotype='{phenotype}',
        p_value_filtering_threshold=p_value_filtering_threshold,
        gff_file=gff_file,
        gff_category=gff_category,
        gff_key=gff_key
    output:
        out_file=os.path.join(os.path.abspath(output_folder),'generate_phenotype_distribution_with_statistical_testing',input_sample + '_{chromosome}__{phenotype}.txt')
    log:
        os.path.join(os.path.abspath(output_folder),'generate_phenotype_distribution_with_statistical_testing_log',input_sample + '_{chromosome}__{phenotype}.log')
    threads: threads
    resources:
        memory=memory
    shell:
        """
        python3 {workflow_path}/scripts/python/generate_phenotype_distribution_with_statistical_testing.py \
        -i {input.in_file} \
        -c {input.accession_mapping_file} \
        -ak {params.accession_key} \
        -pk {params.phenotype_key} \
        -r {input.phenotype_file} \
        -ph {params.phenotype} \
        -pt {params.p_value_filtering_threshold} \
        -g {params.gff_file} \
        -gc "{params.gff_category}" \
        -gk "{params.gff_key}" \
        -o {output.out_file} 2>&1 > {log}
        """
