# PhenoDist

<!-- badges: start -->
<!-- badges: end -->

The PhenoDist pipeline is built for generating phenotype distributions for alleles in variant positions and utilizing statistical methods to test for variant position significance.

## Requirements

In order to run the PhenoDist pipeline, users need to install Miniconda and prepare the Miniconda environment in their computing systems.

The required software, programming languages, and packages include:

```
bwa>=0.7.17
gatk4>=4.4.0.0
samtools>=1.6
htslib>=1.3
python>=3.12
snakemake>=8.4
scipy>=1.12
numpy>=1.26
pandas>=2.2
Beagle>=5.2
SnpEff>=4.3
``` 

Miniconda can be downloaded from [https://docs.anaconda.com/free/miniconda/](https://docs.anaconda.com/free/miniconda/).

For example, if users plan to install Miniconda3 Linux 64-bit, the wget tool can be used to download the Miniconda.

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

To install Miniconda in a server or cluster, users can use the command below.

Please remember to replace the _<installation_shell_script>_ with the actual Miniconda installation shell script. In our case, it is **Miniconda3-latest-Linux-x86_64.sh**.

Please also remember to replace the _<desired_new_directory>_ with an actual directory absolute path.

```
chmod 777 -R <installation_shell_script>
./<installation_shell_script> -b -u -p <desired_new_directory>
rm -rf <installation_shell_script>
```

After installing Miniconda, initialization of Miniconda for bash shell can be done using the command below.

Please also remember to replace the _<desired_new_directory>_ with an actual directory absolute path.

```
<desired_new_directory>/bin/conda init bash
```

Installation of the Miniconda is required, and Miniconda environment needs to be activated every time before running the PhenoDist pipeline.

Write a Conda configuration file (.condarc) before creating a Conda environment:

```
nano ~/.condarc
```

Put the following text into the Conda configuration file (make sure you change _envs_dirs_ and _pkgs_dirs_) then save the file.

Please make sure not use tab in this yaml file, use 4 spaces instead.

Please make sure to replace _/new/path/to/_ with an actual directory absolute path.

```
envs_dirs:
    - /new/path/to/miniconda/envs
pkgs_dirs:
    - /new/path/to/miniconda/pkgs
channels:
    - conda-forge
    - bioconda
    - defaults
```

Create a Conda environment by specifying all required packages (option 1).

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda create -n <conda_environment_name> bioconda::gatk4 bioconda::samtools bioconda::bcftools bioconda::htslib \
bioconda::bedtools bioconda::bwa bioconda::snakemake bioconda::snakemake-executor-plugin-cluster-generic \
conda-forge::numpy conda-forge::pandas
```

Create a Conda environment by using a yaml environment file (option 2).

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda create --name <conda_environment_name> --file PhenoDist-environment.yml
```

Create a Conda environment by using a explicit specification file (option 3).

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda create --name <conda_environment_name> --file PhenoDist-spec-file.txt
```

Activate Conda environment using conda activate command. 

This step is required every time before running PhenoDist pipeline.

Please make sure to replace the _<conda_environment_name>_ with an environment name of your choice.

```
conda activate <conda_environment_name>
```

## Installation

You can install the PhenoDist from [Github](https://github.com/yenon118/PhenoDist.git) with:

```
git clone https://github.com/yenon118/PhenoDist.git
```

Create a "tools" folder within the PhenoDist pipeline clone. 

```
cd PhenoDist
mkdir -p tools
cd tools
```

The Beagle imputation tool can be downloaded from [https://faculty.washington.edu/browning/beagle/beagle.html](https://faculty.washington.edu/browning/beagle/beagle.html) and placed inside the "tools" folder.

The SnpEff functional effect prediction tool can be downloaded from [https://pcingola.github.io/SnpEff/download/](https://pcingola.github.io/SnpEff/download/) and placed in the "tools" folder. 

## Usage

#### Write a configuration file in json format

Please save the file with .json extension.

```
{
	"project_name": "Test_Soy1066",
	"workflow_path": "/storage/htc/joshilab/yenc/projects/PhenoDist",
	"input_files": [
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr01.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr02.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr03.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr04.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr05.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr06.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr07.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr08.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr09.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr10.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr11.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr12.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr13.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr14.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr15.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr16.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr17.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr18.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr19.vcf.gz",
		"/storage/htc/joshilab/yenc/projects/PhenoDist/data/test_Soy1066/unimpute_data/test_Soy1066_Chr20.vcf.gz"
	],
	"chromosomes": [
		"Chr01",
		"Chr02",
		"Chr03",
		"Chr04",
		"Chr05",
		"Chr06",
		"Chr07",
		"Chr08",
		"Chr09",
		"Chr10",
		"Chr11",
		"Chr12",
		"Chr13",
		"Chr14",
		"Chr15",
		"Chr16",
		"Chr17",
		"Chr18",
		"Chr19",
		"Chr20"
	],
	"beagle_window": 40,
	"genome_version": "Wm82.a2.v1",
	"accession_mapping_file": "/storage/htc/joshilab/yenc/projects/PhenoDist/data/act_Soy2939_Accession_Mapping.tsv",
	"accession_key": "Accession",
	"phenotype_key": "GRIN_Accession",
	"phenotype_file": "/storage/htc/joshilab/yenc/projects/PhenoDist/data/act_Soy2939_Phenotype_Data.tsv",
	"phenotypes": [
		"HEIGHT",
		"PODCOLOR"
	],
	"output_folder": "/storage/htc/joshilab/yenc/projects/PhenoDist/output/Test_Soy1066/",
	"memory": 100,
	"threads": 4
}
```

#### Run workflow with the Snakemake workflow management system

```
snakemake -j NUMBER_OF_JOBS --configfile CONFIGURATION_FILE --snakefile SNAKEMAKE_FILE

Mandatory Positional Argumants:
    NUMBER_OF_JOBS                          - the number of jobs
    CONFIGURATION_FILE                      - a configuration file
    SNAKEMAKE_FILE                          - the PhenoDist.smk file that sit inside this repository 
```

## Examples

Below are some fundamental examples illustrating the usage of the PhenoDist pipeline.

Please adjust _/path/to/_ to an actual directory absolute path.

**Examples of running without an executor.**

```
cd /path/to/PhenoDist

snakemake -pj 3 --configfile lewis_slurm_inputs.json --snakefile PhenoDist.smk
```

**Examples of running with an executor.**

Snakemake version >= 8.0.0.

```
cd /path/to/PhenoDist

snakemake --executor cluster-generic \
--cluster-generic-submit-cmd "sbatch --account=xulab --time=0-02:00 \
--partition=Lewis,BioCompute,hpc5,General --mem=64G" \
--jobs 25 --latency-wait 60 \
--configfile lewis_slurm_inputs.json \
--snakefile PhenoDist.smk
```

Snakemake version < 8.0.0.

```
cd /path/to/PhenoDist

snakemake --cluster "sbatch --account=xulab --time=0-02:00 \
--partition=Lewis,BioCompute,hpc5,General --mem=64G" \
--jobs 25 --latency-wait 60 \
--configfile lewis_slurm_inputs.json \
--snakefile PhenoDist.smk
```
