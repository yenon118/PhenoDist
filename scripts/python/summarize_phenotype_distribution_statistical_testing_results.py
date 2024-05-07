#!/usr/bin/env python3

import sys
import os
import re
import argparse
import pathlib
import gzip


def main(args):
    #######################################################################
    # Get arguments
    #######################################################################
    input_file_paths = args.input_file
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
    # Write output file header
    #######################################################################
    with open(output_file_path, 'w') as writer:
        writer.write(
            "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                'Gene',
                'Phenotype',
                'Phenotype_Data_Type',
                'Test_Method',
                'Minimum_Test_P_Value',
                'Maximum_Test_P_Value'
            )
        )

    #######################################################################
    # Read input file
    #######################################################################
    output_dict = {}

    for input_file_path in input_file_paths:
        with open(input_file_path, "r") as reader:
            header = reader.readline()
            for line in reader:
                line_array = str(line).strip("\n").strip("\r").strip("\r\n").split("\t")

                gene = line_array[2]
                phenotype = line_array[5]
                phenotype_data_type = line_array[6]
                test_method = line_array[7]
                res_pvalue = line_array[14]

                if gene not in output_dict:
                    output_dict[gene] = {}

                if phenotype not in output_dict[gene]:
                    output_dict[gene][phenotype] = {}

                if test_method not in output_dict[gene][phenotype]:
                    output_dict[gene][phenotype][test_method] = {
                        'Phenotype': phenotype,
                        'Phenotype_Data_Type': phenotype_data_type,
                        'Test_Method': test_method,
                        'Minimum_Test_P_Value': '',
                        'Maximum_Test_P_Value': ''
                    }

                if res_pvalue != '':
                    if output_dict[gene][phenotype][test_method]['Minimum_Test_P_Value'] == '':
                        output_dict[gene][phenotype][test_method]['Minimum_Test_P_Value'] = res_pvalue
                    else:
                        if float(res_pvalue) < float(output_dict[gene][phenotype][test_method]['Minimum_Test_P_Value']):
                            output_dict[gene][phenotype][test_method]['Minimum_Test_P_Value'] = res_pvalue
                    if output_dict[gene][phenotype][test_method]['Maximum_Test_P_Value'] == '':
                        output_dict[gene][phenotype][test_method]['Maximum_Test_P_Value'] = res_pvalue
                    else:
                        if float(res_pvalue) > float(output_dict[gene][phenotype][test_method]['Maximum_Test_P_Value']):
                            output_dict[gene][phenotype][test_method]['Maximum_Test_P_Value'] = res_pvalue

    # Check and write data
    if output_dict:
        with open(output_file_path, 'a') as writer:
            for gene in output_dict:
                for phenotype in output_dict[gene]:
                    for test_method in output_dict[gene][phenotype]:
                        writer.write(
                            "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                gene,
                                phenotype,
                                output_dict[gene][phenotype][test_method]['Phenotype_Data_Type'],
                                test_method,
                                output_dict[gene][phenotype][test_method]['Minimum_Test_P_Value'],
                                output_dict[gene][phenotype][test_method]['Maximum_Test_P_Value']
                            )
                        )


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='summarize_phenotype_distribution_statistical_testing_results', description='summarize_phenotype_distribution_statistical_testing_results')

    parser.add_argument('-i', '--input_file', help='Input file', action='append', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
