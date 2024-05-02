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

    for input_file_path in input_file_paths:
        with open(input_file_path, "r") as reader:
            header = reader.readline()
            for line in reader:
                output_array.append(str(line))
                # Check and write data
                if (len(output_array) > chunksize):
                    with open(output_file_path, 'a') as writer:
                        writer.write("".join(output_array))
                        output_array.clear()

    # Check and write data
    if (len(output_array) > 0):
        with open(output_file_path, 'a') as writer:
            writer.write("".join(output_array))
            output_array.clear()


if __name__ == "__main__":
    #######################################################################
    # Parse arguments
    #######################################################################
    parser = argparse.ArgumentParser(prog='combine_phenotype_distribution_statistical_testing_results', description='combine_phenotype_distribution_statistical_testing_results')

    parser.add_argument('-i', '--input_file', help='Input file', action='append', type=pathlib.Path, required=True)
    parser.add_argument('-o', '--output_file', help='Output file', type=pathlib.Path, required=True)

    args = parser.parse_args()

    #######################################################################
    # Call main function
    #######################################################################
    main(args)
