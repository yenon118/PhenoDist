#!/usr/bin/env python3

import sys
import os
import re
import math
import argparse
import pathlib
import gzip
import pprint
import scipy
import warnings

from itertools import combinations, permutations

import numpy as np
import pandas as pd

warnings.filterwarnings("error")


# Process line in VCF file
def process_line(header_array, line_array, gff_dict, accession_mapping_dict, accession_key, phenotype_key, phenotype_dict, phenotype, p_value_filtering_threshold, output_array):
    chromosome = line_array[0]
    position = line_array[1]
    reference_allele = line_array[3]
    alternate_alleles = line_array[4]

    alleles = re.split(',', reference_allele) + re.split(',', alternate_alleles)

    accession_array = header_array[9:]

    genotype_index_array = [int(str(re.split('/|\\|', genotype_index)[0]).strip()) for genotype_index in line_array[9:]]

    genotype_array = [alleles[genotype_index] for genotype_index in genotype_index_array]

    category_array = ['Ref' if genotype_index == 0 else 'Alt' for genotype_index in genotype_index_array]

    category_index_array = [0 if genotype_index == 0 else 1 for genotype_index in genotype_index_array]

    gene_array = []
    if chromosome in gff_dict:
        for gene in gff_dict[chromosome]:
            if int(gff_dict[chromosome][gene]['Start']) <= int(position) <= int(gff_dict[chromosome][gene]['End']):
                gene_array.append(gene)

    if phenotype in phenotype_dict:
        if phenotype_dict[phenotype]['DataType'] != '':
            phenotype_data_type = phenotype_dict[phenotype]['DataType']

            # Generate genotype_phenotype_dict: genotype_phenotype_dict[`allele`] = `alleles and phenotype data`
            genotype_phenotype_dict = {}
            for i in range(len(genotype_array)):
                if genotype_array[i] not in genotype_phenotype_dict.keys():
                    genotype_phenotype_dict[genotype_array[i]] = []
                if accession_array[i] in accession_mapping_dict.keys():
                    if phenotype in phenotype_dict.keys():
                        if accession_mapping_dict[accession_array[i]] in phenotype_dict[phenotype]['Data'].keys():
                            if phenotype_dict[phenotype]['Data'][accession_mapping_dict[accession_array[i]]] != '':
                                new_dict = {
                                    'Accession': accession_array[i],
                                    'Accession_Mapping': accession_mapping_dict[accession_array[i]],
                                    'Phenotype_Measurement': phenotype_dict[phenotype]['Data'][accession_mapping_dict[accession_array[i]]]
                                }
                                genotype_phenotype_dict[genotype_array[i]].append(new_dict)

            # Create allele_phenotype_dict for alleles and phenotype values
            allele_phenotype_dict = {}
            for allele in genotype_phenotype_dict.keys():
                allele_phenotype_dict[allele] = []
                for i in range(len(genotype_phenotype_dict[allele])):
                    allele_phenotype_dict[allele] = allele_phenotype_dict[allele] + genotype_phenotype_dict[allele][i]['Phenotype_Measurement']

            if len(allele_phenotype_dict.keys()) > 1:
                # Generate allele combinations and/or permutations
                allele_combination_list = list(combinations(allele_phenotype_dict.keys(), 2))
                # allele_permutation_list = list(permutations(allele_phenotype_dict.keys(), 2))

                if phenotype_data_type == 'float':
                    for allele_combination in allele_combination_list:

                        allele_combination_1 = allele_combination[0]
                        allele_combination_2 = allele_combination[1]
                        test_method = ''
                        phenotype_category_1 = ''
                        phenotype_category_2 = ''
                        res_shapiro = None
                        res_shapiro_statistic = ''
                        res_shapiro_pvalue = ''
                        res_mannwhitneyu = None
                        res_mannwhitneyu_statistic = ''
                        res_mannwhitneyu_pvalue = ''
                        mannwhitneyu_accession_count_1 = ''
                        mannwhitneyu_accession_count_2 = ''
                        res_ttest = None
                        res_ttest_statistic = ''
                        res_ttest_pvalue = ''
                        ttest_accession_count_1 = ''
                        ttest_accession_count_2 = ''
                        res_chi2_contingency = None
                        res_chi2_contingency_statistic = ''
                        res_chi2_contingency_pvalue = ''
                        chi2_contingency_accession_count_1 = ''
                        chi2_contingency_accession_count_2 = ''
                        res_fisher_exact_array = []
                        res_fisher_exact_statistic = ''
                        res_fisher_exact_pvalue = ''
                        fisher_exact_accession_count = ''

                        # Make sure data is valid for shapiro normality testing
                        if allele_phenotype_dict[allele_combination[0]] and allele_phenotype_dict[allele_combination[1]]:
                            if len(allele_phenotype_dict[allele_combination[0]] + allele_phenotype_dict[allele_combination[1]]) > 2:
                                # Run shapiro test for normality testing
                                try:
                                    res_shapiro = scipy.stats.shapiro(allele_phenotype_dict[allele_combination[0]] + allele_phenotype_dict[allele_combination[1]])
                                    if math.isfinite(res_shapiro.statistic):
                                        res_shapiro_statistic = res_shapiro.statistic
                                    else:
                                        res_shapiro_statistic = ''
                                    if math.isfinite(res_shapiro.pvalue):
                                        res_shapiro_pvalue = res_shapiro.pvalue
                                    else:
                                        res_shapiro_pvalue = ''
                                except Exception as e:
                                    res_shapiro = None
                                if res_shapiro is not None:
                                    if res_shapiro_pvalue != '':
                                        if float(res_shapiro_pvalue) < float(p_value_filtering_threshold):
                                            # Not normally distributed, run Mann-Whitney U rank test
                                            try:
                                                res_mannwhitneyu = scipy.stats.mannwhitneyu(
                                                    allele_phenotype_dict[allele_combination[0]],
                                                    allele_phenotype_dict[allele_combination[1]]
                                                )
                                                if math.isfinite(res_mannwhitneyu.statistic):
                                                    res_mannwhitneyu_statistic = res_mannwhitneyu.statistic
                                                else:
                                                    res_mannwhitneyu_statistic = ''
                                                if math.isfinite(res_mannwhitneyu.pvalue):
                                                    res_mannwhitneyu_pvalue = res_mannwhitneyu.pvalue
                                                else:
                                                    res_mannwhitneyu_pvalue = ''
                                                # Specify test method
                                                test_method = 'Mann-Whitney U Rank Test'
                                                # Get accession counts
                                                mannwhitneyu_accession_count_1 = len(allele_phenotype_dict[allele_combination[0]])
                                                mannwhitneyu_accession_count_2 = len(allele_phenotype_dict[allele_combination[1]])
                                            except RuntimeWarning as e:
                                                res_mannwhitneyu = None
                                            except Exception as e:
                                                res_mannwhitneyu = None

                                            if gene_array:
                                                if len(gene_array) > 0:
                                                    for gene in gene_array:
                                                        if res_mannwhitneyu:
                                                            if res_mannwhitneyu_pvalue != '':
                                                                if float(res_mannwhitneyu_pvalue) < float(p_value_filtering_threshold):
                                                                    output_array.append(
                                                                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                                            chromosome,
                                                                            position,
                                                                            gene,
                                                                            allele_combination_1,
                                                                            allele_combination_2,
                                                                            phenotype,
                                                                            phenotype_data_type,
                                                                            test_method,
                                                                            phenotype_category_1,
                                                                            phenotype_category_2,
                                                                            res_shapiro_statistic,
                                                                            res_shapiro_pvalue,
                                                                            res_mannwhitneyu_statistic,
                                                                            res_mannwhitneyu_pvalue,
                                                                            mannwhitneyu_accession_count_1,
                                                                            mannwhitneyu_accession_count_2,
                                                                            res_ttest_statistic,
                                                                            res_ttest_pvalue,
                                                                            ttest_accession_count_1,
                                                                            ttest_accession_count_2,
                                                                            res_chi2_contingency_statistic,
                                                                            res_chi2_contingency_pvalue,
                                                                            chi2_contingency_accession_count_1,
                                                                            chi2_contingency_accession_count_2,
                                                                            res_fisher_exact_statistic,
                                                                            res_fisher_exact_pvalue,
                                                                            fisher_exact_accession_count
                                                                        )
                                                                    )
                                        else:
                                            # Normally distributed, run t-test
                                            try:
                                                res_ttest = scipy.stats.ttest_ind(
                                                    allele_phenotype_dict[allele_combination[0]],
                                                    allele_phenotype_dict[allele_combination[1]]
                                                )
                                                if math.isfinite(res_ttest.statistic):
                                                    res_ttest_statistic = res_ttest.statistic
                                                else:
                                                    res_ttest_statistic = ''
                                                if math.isfinite(res_ttest.pvalue):
                                                    res_ttest_pvalue = res_ttest.pvalue
                                                else:
                                                    res_ttest_pvalue = ''
                                                # Specify test method
                                                test_method = 'T-test'
                                                # Get accession counts
                                                ttest_accession_count_1 = len(allele_phenotype_dict[allele_combination[0]])
                                                ttest_accession_count_2 = len(allele_phenotype_dict[allele_combination[1]])
                                            except RuntimeWarning as e:
                                                res_ttest = None
                                            except Exception as e:
                                                res_ttest = None

                                            if gene_array:
                                                if len(gene_array) > 0:
                                                    for gene in gene_array:
                                                        if res_ttest:
                                                            if res_ttest_pvalue != '':
                                                                if float(res_ttest_pvalue) < float(p_value_filtering_threshold):
                                                                    output_array.append(
                                                                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                                            chromosome,
                                                                            position,
                                                                            gene,
                                                                            allele_combination_1,
                                                                            allele_combination_2,
                                                                            phenotype,
                                                                            phenotype_data_type,
                                                                            test_method,
                                                                            phenotype_category_1,
                                                                            phenotype_category_2,
                                                                            res_shapiro_statistic,
                                                                            res_shapiro_pvalue,
                                                                            res_mannwhitneyu_statistic,
                                                                            res_mannwhitneyu_pvalue,
                                                                            mannwhitneyu_accession_count_1,
                                                                            mannwhitneyu_accession_count_2,
                                                                            res_ttest_statistic,
                                                                            res_ttest_pvalue,
                                                                            ttest_accession_count_1,
                                                                            ttest_accession_count_2,
                                                                            res_chi2_contingency_statistic,
                                                                            res_chi2_contingency_pvalue,
                                                                            chi2_contingency_accession_count_1,
                                                                            chi2_contingency_accession_count_2,
                                                                            res_fisher_exact_statistic,
                                                                            res_fisher_exact_pvalue,
                                                                            fisher_exact_accession_count
                                                                        )
                                                                    )
                elif phenotype_data_type == 'string':
                    for allele_combination in allele_combination_list:

                        allele_combination_1 = allele_combination[0]
                        allele_combination_2 = allele_combination[1]
                        test_method = ''
                        phenotype_category_1 = ''
                        phenotype_category_2 = ''
                        res_shapiro = None
                        res_shapiro_statistic = ''
                        res_shapiro_pvalue = ''
                        res_mannwhitneyu = None
                        res_mannwhitneyu_statistic = ''
                        res_mannwhitneyu_pvalue = ''
                        mannwhitneyu_accession_count_1 = ''
                        mannwhitneyu_accession_count_2 = ''
                        res_ttest = None
                        res_ttest_statistic = ''
                        res_ttest_pvalue = ''
                        ttest_accession_count_1 = ''
                        ttest_accession_count_2 = ''
                        res_chi2_contingency = None
                        res_chi2_contingency_statistic = ''
                        res_chi2_contingency_pvalue = ''
                        chi2_contingency_accession_count_1 = ''
                        chi2_contingency_accession_count_2 = ''
                        res_fisher_exact_array = []
                        res_fisher_exact_statistic = ''
                        res_fisher_exact_pvalue = ''
                        fisher_exact_accession_count = ''

                        # Grab all distinct phenotype categories
                        unique_phenotype_category_array = []
                        for i in range(len(allele_combination)):
                            for phenotype_category in allele_phenotype_dict[allele_combination[i]]:
                                if phenotype_category not in unique_phenotype_category_array:
                                    if phenotype_category != '':
                                        unique_phenotype_category_array.append(phenotype_category)

                        if unique_phenotype_category_array:
                            if len(unique_phenotype_category_array) > 1:
                                # Generate allele_phenotype_category_count_dict
                                allele_phenotype_category_count_dict = {}
                                for i in range(len(allele_combination)):
                                    allele_phenotype_category_count_dict[allele_combination[i]] = {}
                                    for phenotype_category in unique_phenotype_category_array:
                                        allele_phenotype_category_count_dict[allele_combination[i]][phenotype_category] = 0

                                # Get count for allele and phenotype category
                                for i in range(len(allele_combination)):
                                    if allele_combination[i] in allele_phenotype_dict:
                                        if allele_phenotype_dict[allele_combination[i]]:
                                            for phenotype_category in allele_phenotype_dict[allele_combination[i]]:
                                                allele_phenotype_category_count_dict[allele_combination[i]][phenotype_category] += 1

                                # Run chi-square test
                                try:
                                    res_chi2_contingency = scipy.stats.chi2_contingency([
                                        list(allele_phenotype_category_count_dict[allele_combination[0]].values()),
                                        list(allele_phenotype_category_count_dict[allele_combination[1]].values())
                                    ])
                                    if math.isfinite(res_chi2_contingency.statistic):
                                        res_chi2_contingency_statistic = res_chi2_contingency.statistic
                                    else:
                                        res_chi2_contingency_statistic = ''
                                    if math.isfinite(res_chi2_contingency.pvalue):
                                        res_chi2_contingency_pvalue = res_chi2_contingency.pvalue
                                    else:
                                        res_chi2_contingency_pvalue = ''
                                    # Specify test method
                                    test_method = 'Chi-square Test'
                                    # Get accession counts
                                    chi2_contingency_accession_count_1 = ','.join([str(i) for i in list(allele_phenotype_category_count_dict[allele_combination[0]].values())])
                                    chi2_contingency_accession_count_2 = ','.join([str(i) for i in list(allele_phenotype_category_count_dict[allele_combination[1]].values())])
                                except Exception as e:
                                    res_chi2_contingency = None

                                if gene_array:
                                    if len(gene_array) > 0:
                                        for gene in gene_array:
                                            if res_chi2_contingency:
                                                if res_chi2_contingency_pvalue != '':
                                                    if float(res_chi2_contingency_pvalue) < float(p_value_filtering_threshold):
                                                        output_array.append(
                                                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                                chromosome,
                                                                position,
                                                                gene,
                                                                allele_combination_1,
                                                                allele_combination_2,
                                                                phenotype,
                                                                phenotype_data_type,
                                                                test_method,
                                                                phenotype_category_1,
                                                                phenotype_category_2,
                                                                res_shapiro_statistic,
                                                                res_shapiro_pvalue,
                                                                res_mannwhitneyu_statistic,
                                                                res_mannwhitneyu_pvalue,
                                                                mannwhitneyu_accession_count_1,
                                                                mannwhitneyu_accession_count_2,
                                                                res_ttest_statistic,
                                                                res_ttest_pvalue,
                                                                ttest_accession_count_1,
                                                                ttest_accession_count_2,
                                                                res_chi2_contingency_statistic,
                                                                res_chi2_contingency_pvalue,
                                                                chi2_contingency_accession_count_1,
                                                                chi2_contingency_accession_count_2,
                                                                res_fisher_exact_statistic,
                                                                res_fisher_exact_pvalue,
                                                                fisher_exact_accession_count
                                                            )
                                                        )

                                # Generate phenotype category combinations and/or permutations
                                phenotype_category_combination_list = list(combinations(unique_phenotype_category_array, 2))

                                # Run fisher exact test
                                for phenotype_category_combination in phenotype_category_combination_list:
                                    try:
                                        res_fisher_exact = scipy.stats.fisher_exact(
                                            [
                                                [
                                                    allele_phenotype_category_count_dict[allele_combination[0]][phenotype_category_combination[0]],
                                                    allele_phenotype_category_count_dict[allele_combination[0]][phenotype_category_combination[1]]
                                                ],
                                                [
                                                    allele_phenotype_category_count_dict[allele_combination[1]][phenotype_category_combination[0]],
                                                    allele_phenotype_category_count_dict[allele_combination[1]][phenotype_category_combination[1]]
                                                ]
                                            ]
                                        )
                                        if math.isfinite(res_fisher_exact.statistic):
                                            res_fisher_exact_statistic = res_fisher_exact.statistic
                                        else:
                                            res_fisher_exact_statistic = ''
                                        if math.isfinite(res_fisher_exact.pvalue):
                                            res_fisher_exact_pvalue = res_fisher_exact.pvalue
                                        else:
                                            res_fisher_exact_pvalue = ''
                                    except Exception as e:
                                        res_fisher_exact = None

                                    if res_fisher_exact is not None:
                                        if res_fisher_exact_pvalue != '':
                                            if float(res_fisher_exact_pvalue) < float(p_value_filtering_threshold):
                                                res_fisher_exact_array.append({
                                                    'Allele_1': allele_combination[0],
                                                    'Allele_2': allele_combination[1],
                                                    'Test_Method': 'Fisher Exact Test',
                                                    'Phenotype_Category_1': phenotype_category_combination[0],
                                                    'Phenotype_Category_2': phenotype_category_combination[1],
                                                    'Chi_Square_Test_Statistics': '',
                                                    'Chi_Square_Test_P_Value': '',
                                                    'Chi_Square_Test_Accession_Count_1': '',
                                                    'Chi_Square_Test_Accession_Count_2': '',
                                                    'Fisher_Exact_Test_Statistics': res_fisher_exact_statistic,
                                                    'Fisher_Exact_Test_P_Value': res_fisher_exact_pvalue,
                                                    'Fisher_Exact_Test_Accession_Count': str(allele_phenotype_category_count_dict[allele_combination[0]][
                                                            phenotype_category_combination[0]]) + "," + str(
                                                        allele_phenotype_category_count_dict[allele_combination[0]][
                                                            phenotype_category_combination[1]]) + "," + str(
                                                        allele_phenotype_category_count_dict[allele_combination[1]][
                                                            phenotype_category_combination[0]]) + "," + str(
                                                        allele_phenotype_category_count_dict[allele_combination[1]][
                                                            phenotype_category_combination[1]])
                                                })

                                if gene_array:
                                    if len(gene_array) > 0:
                                        for gene in gene_array:
                                            if res_fisher_exact_array:
                                                if len(res_fisher_exact_array) > 0:
                                                    for i in range(len(res_fisher_exact_array)):
                                                        output_array.append(
                                                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                                                chromosome,
                                                                position,
                                                                gene,
                                                                res_fisher_exact_array[i]['Allele_1'],
                                                                res_fisher_exact_array[i]['Allele_2'],
                                                                phenotype,
                                                                phenotype_data_type,
                                                                res_fisher_exact_array[i]['Test_Method'],
                                                                res_fisher_exact_array[i]['Phenotype_Category_1'],
                                                                res_fisher_exact_array[i]['Phenotype_Category_2'],
                                                                res_shapiro_statistic,
                                                                res_shapiro_pvalue,
                                                                res_mannwhitneyu_statistic,
                                                                res_mannwhitneyu_pvalue,
                                                                mannwhitneyu_accession_count_1,
                                                                mannwhitneyu_accession_count_2,
                                                                res_ttest_statistic,
                                                                res_ttest_pvalue,
                                                                ttest_accession_count_1,
                                                                ttest_accession_count_2,
                                                                res_fisher_exact_array[i]['Chi_Square_Test_Statistics'],
                                                                res_fisher_exact_array[i]['Chi_Square_Test_P_Value'],
                                                                res_fisher_exact_array[i]['Chi_Square_Test_Accession_Count_1'],
                                                                res_fisher_exact_array[i]['Chi_Square_Test_Accession_Count_2'],
                                                                res_fisher_exact_array[i]['Fisher_Exact_Test_Statistics'],
                                                                res_fisher_exact_array[i]['Fisher_Exact_Test_P_Value'],
                                                                res_fisher_exact_array[i]['Fisher_Exact_Test_Accession_Count']
                                                            )
                                                        )


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_path = args.input_file
    accession_mapping_file_path = args.accession_mapping_file
    accession_key = args.accession_key
    phenotype_key = args.phenotype_key
    phenotype_file_path = args.phenotype_file
    phenotype = args.phenotype
    p_value_filtering_threshold = args.p_value_filtering_threshold
    gff_file_path = args.gff_file
    gff_category = args.gff_category
    gff_key = args.gff_key
    output_file_path = args.output_file

    #######################################################################
    # Check if output parent folder exists
    # If not, create the output parent folder
    #######################################################################
    if not output_file_path.parent.exists():
        try:
            output_file_path.parent.mkdir(parents=True)
        except FileNotFoundError as e:
            pass
        except FileExistsError as e:
            pass
        except Exception as e:
            pass
        if not output_file_path.parent.exists():
            sys.exit(1)

    #######################################################################
    # Load GFF file
    #######################################################################
    gff_dict = {}

    with open(gff_file_path, 'r') as reader:
        for line in reader:
            if not line.startswith('#'):
                line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")

                if str(str(line_array[2]).strip()) == gff_category:

                    chromosome = str(line_array[0]).strip()
                    start = str(line_array[3]).strip()
                    end = str(line_array[4]).strip()
                    attributes = str(line_array[8]).strip()

                    gene_id = re.sub('(;.*)', '', re.sub(re.compile(".*" + str(gff_key)), '', attributes))

                    if chromosome not in gff_dict.keys():
                        gff_dict[chromosome] = {}

                    if gene_id not in gff_dict[chromosome].keys():
                        gff_dict[chromosome][gene_id] = {
                            "Start": start,
                            "End": end
                        }
                    else:
                        gff_dict[chromosome][gene_id]['Start'] = start
                        gff_dict[chromosome][gene_id]['End'] = end

    #######################################################################
    # Read accession mapping
    #######################################################################
    accession_mapping_dict = {}

    with open(accession_mapping_file_path, "r") as reader:
        header = reader.readline()
        header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
        header_array = [element.strip("\"") for element in header_array]
        # Find accession key index
        try:
            accession_key_index = header_array.index(accession_key)
        except Exception as e:
            accession_key_index = None
            print("Accession key index cannot be found!!!")
            print(e)
        # Find phenotype key index
        try:
            phenotype_key_index = header_array.index(phenotype_key)
        except Exception as e:
            phenotype_key_index = None
            print("Phenotype key index cannot be found!!!")
            print(e)
        for line in reader:
            line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
            line_array = [element.strip("\"") for element in line_array]
            # Populate accession_mapping_dict: accession_mapping_dict[`accession key`] = `phenotype key`
            if (accession_key_index is not None) and (phenotype_key_index is not None):
                if ((line_array[accession_key_index] is not None) and (line_array[phenotype_key_index] is not None)):
                    if (line_array[accession_key_index] != '') and (line_array[phenotype_key_index] != ''):
                        if line_array[accession_key_index] not in accession_mapping_dict.keys():
                            accession_mapping_dict[line_array[accession_key_index]] = line_array[phenotype_key_index]

    #######################################################################
    # Read phenotype data
    #######################################################################
    phenotype_dict = {}

    with open(phenotype_file_path, "r") as reader:
        header = reader.readline()
        header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
        header_array = [element.strip("\"") for element in header_array]
        # Add column names to phenotype_dict and prepare for data type and data
        for i in range(1, len(header_array)):
            if header_array[i] is not None:
                if header_array[i] != '':
                    if header_array[i] not in phenotype_dict.keys():
                        phenotype_dict[header_array[i]] = {
                            'DataType': '',
                            'Data': {}
                        }
        for line in reader:
            line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
            line_array = [element.strip("\"") for element in line_array]
            # Add accessions and values into phenotype_dict[`phenotype`]['Data']
            for i in range(1, len(line_array)):
                if (line_array[0] is not None) and (line_array[i] is not None):
                    if (line_array[0] != '') and (line_array[i] != ''):
                        if header_array[i] in phenotype_dict.keys():
                            phenotype_value_array = []
                            for element in re.split(',', line_array[i]):
                                if element != '':
                                    phenotype_value_array.append(element)
                            phenotype_dict[header_array[i]]['Data'][line_array[0]] = phenotype_value_array

    #######################################################################
    # Check and convert phenotype data
    #######################################################################

    for pheno in phenotype_dict.keys():
        if phenotype_dict[pheno]['Data']:
            is_float_flag = True
            # Check if data is in float data type
            for acc in phenotype_dict[pheno]['Data'].keys():
                try:
                    converted_value = [float(element) for element in phenotype_dict[pheno]['Data'][acc]]
                except Exception as e:
                    is_float_flag = False
                    pass
            # Convert all values to float data type if data type is float else keep as string data type
            if is_float_flag:
                phenotype_dict[pheno]['DataType'] = 'float'
                for acc in phenotype_dict[pheno]['Data'].keys():
                    phenotype_dict[pheno]['Data'][acc] = [float(element) for element in phenotype_dict[pheno]['Data'][acc]]
            else:
                phenotype_dict[pheno]['DataType'] = 'string'

    #######################################################################
    # Write output file header
    #######################################################################
    with open(output_file_path, 'w') as writer:
        writer.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                'Chromosome',
                'Position',
                'Gene',
                'Allele_1',
                'Allele_2',
                'Phenotype',
                'Phenotype_Data_Type',
                'Test_Method',
                'Phenotype_Category_1',
                'Phenotype_Category_2',
                'Shapiro_Test_Statistics',
                'Shapiro_Test_P_Value',
                'Mann_Whitney_U_Rank_Test_Statistics',
                'Mann_Whitney_U_Rank_Test_P_Value',
                'Mann_Whitney_U_Rank_Test_Accession_Count_1',
                'Mann_Whitney_U_Rank_Test_Accession_Count_2',
                'T_Test_Statistics',
                'T_Test_P_Value',
                'T_Test_Accession_Count_1',
                'T_Test_Accession_Count_2',
                'Chi_Square_Test_Statistics',
                'Chi_Square_Test_P_Value',
                'Chi_Square_Test_Accession_Count_1',
                'Chi_Square_Test_Accession_Count_2',
                'Fisher_Exact_Test_Statistics',
                'Fisher_Exact_Test_P_Value',
                'Fisher_Exact_Test_Accession_Count'
            )
        )

    #######################################################################
    # Read input file
    #######################################################################
    chunksize = 10000
    output_array = []

    if str(input_file_path).endswith('gz'):
        with gzip.open(input_file_path, 'rt') as reader:
            header = ""
            while not header.strip().startswith("#CHROM"):
                header = reader.readline()
                header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
            for line in reader:
                line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
                process_line(
                    header_array,
                    line_array,
                    gff_dict,
                    accession_mapping_dict,
                    accession_key,
                    phenotype_key,
                    phenotype_dict,
                    phenotype,
                    p_value_filtering_threshold,
                    output_array
                )
                # Check and write data
                if len(output_array) > chunksize:
                    with open(output_file_path, 'a') as writer:
                        writer.write("".join(output_array))
                        output_array.clear()
    else:
        with open(input_file_path, "r") as reader:
            header = ""
            while not header.strip().startswith("#CHROM"):
                header = reader.readline()
                header_array = str(header).strip("\n").strip("\r").strip("\r\n").split("\t")
            for line in reader:
                line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")
                process_line(
                    header_array,
                    line_array,
                    gff_dict,
                    accession_mapping_dict,
                    accession_key,
                    phenotype_key,
                    phenotype_dict,
                    phenotype,
                    p_value_filtering_threshold,
                    output_array
                )
                # Check and write data
                if len(output_array) > chunksize:
                    with open(output_file_path, 'a') as writer:
                        writer.write("".join(output_array))
                        output_array.clear()

    # Check and write data
    if len(output_array) > 0:
        with open(output_file_path, 'a') as writer:
            writer.write("".join(output_array))
            output_array.clear()


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='phenotype_to_genotype_statistical_testing', description='phenotype_to_genotype_statistical_testing')

    parser.add_argument('-i', '--input_file', help='Input file', type=pathlib.Path, required=True)
    parser.add_argument('-c', '--accession_mapping_file', help='Accession mapping file', type=pathlib.Path, required=True)
    parser.add_argument('-ak', '--accession_key', help='Accession mapping file key column for accessions in VCF file', type=str, required=True)
    parser.add_argument('-pk', '--phenotype_key', help='Accession mapping file key column for accessions in phenotype file', type=str, required=True)
    parser.add_argument('-r', '--phenotype_file', help='Phenotype file', type=pathlib.Path, required=True)
    parser.add_argument('-ph', '--phenotype', help='Phenotype', type=str, required=True)
    parser.add_argument('-pt', '--p_value_filtering_threshold', help='P-value filtering threshold', type=float, default=0.05)
    parser.add_argument('-g', '--gff_file', help='GFF file', type=pathlib.Path, required=True)
    parser.add_argument('-gc', '--gff_category', help='Gff category', type=str, default='gene')
    parser.add_argument('-gk', '--gff_key', help='Gff key', type=str, default='Name=')
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
