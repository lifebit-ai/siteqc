# lifebit-ai/siteqc: Usage

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run lifebit-ai/siteqc --input ..
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull lifebit-ai/siteqc
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/siteqc releases page](https://github.com/nf-core/siteqc/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Required SiteQC pipeline arguments
### `--input`
File with list of full paths to bcf files and their indexes. Bcf files can be compressed but in a readable for bcftools format.

Example argument:
```
--input s3://lifebit-featured-datasets/projects/gel/siteqc/input.csv
```
Example file content:
```
bcf,index
s3://lifebit-featured-datasets/projects/gel/siteqc/test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz,s3://lifebit-featured-datasets/projects/gel/siteqc/test_all_chunks_merged_norm_chr10_52955340_55447336.bcf.gz.csi
```
### Triodata files
Pipeline need information about participants and ther trio/family relations. This information can be provided in two ways.
1. Default way - providing file with list of samples and file defining the trios. These files will be fed in triodata_define process that will generate keep and fam files needed for the pipeline. Required parameters: `--included_samples`, `--triofile`. 
2. Optional way - provide keep and fam files directly. This will make pipeline skip triodata_define process and use provided files instead. Required parameters: `--triodata_keep_pheno`, `--triodata_fam`.

In a nutshel, either one or another pair of files is required to run the pipeline.
#### `--included_samples`
File with a list of samples to be included in analysis. (Also refered to as participants, participant IDs, platekeys)
Input to triodata_define process.

Example argument:
```
--included_samples s3://lifebit-featured-datasets/projects/gel/siteqc/sampleList.txt
```
Example file content:
```
HG002
sample_1
sample_2
sample_3
HG003
HG004
```
#### `--triofile` 
File describing the family trios of participants in bcf files.
Input to triodata_define process.
Example argument:
```
--triofile s3://lifebit-featured-datasets/projects/gel/siteqc/extended_interpreted_trios_grch38.tsv
```
Example file content:
```
trio_id cohort_id   family_id   interpreation_version   sampleId    participantId   trio    pedigreeId  motherId    fatherId    isProband   affectionStatus sex personKaryotypicSex yearOfBirth lifeStatus  adoptedStatus   consanguineousParents   participantQCState  vcf_id  vcf_size    genome  denovo_filename
trio_1  cohort1 family_1    xxxx    sample_1    pid1    Father  xxxxx   1234    4567    FALSE   UNAFFECTED  MALE    UNKNOWN xxxx    ALIVE   notadopted  unknown passedMedicalReviewReadyForInterpretation   1   3   GRCh38  xxxx
trio_1  cohort1 family_1    xxxx    sample_2    pid2    Mother  xxxxx   2345    5678    FALSE   UNAFFECTED  FEMALE  UNKNOWN xxxx    ALIVE   notadopted  unknown passedMedicalReviewReadyForInterpretation   3   3   GRCh38  xxxx
trio_1  cohort1 family_1    xxxx    sample_3    pid3    Offspring   xxxxx   3456    6789    TRUE    AFFECTED    MALE    UNKNOWN xxxx    ALIVE   notadopted  unknown passedMedicalReviewReadyForInterpretation   2   3   GRCh38  xxxx
trio_2  cohort1 family_2    xxxx    HG002    pid4    Father  xxxxx   1234    4567    FALSE   UNAFFECTED  MALE    UNKNOWN xxxx    ALIVE   notadopted  unknown passedMedicalReviewReadyForInterpretation   1   3   GRCh38  xxxx
trio_2  cohort1 family_2    xxxx    HG003    pid5    Mother  xxxxx   2345    5678    FALSE   UNAFFECTED  FEMALE  UNKNOWN xxxx    ALIVE   notadopted  unknown passedMedicalReviewReadyForInterpretation   3   3   GRCh38  xxxx
trio_2  cohort1 family_2    xxxx    HG004    pid6    Offspring   xxxxx   3456    6789    TRUE    AFFECTED    MALE    UNKNOWN xxxx    ALIVE   notadopted  unknown passedMedicalReviewReadyForInterpretation   2   3   GRCh38  xxxx
```
#### `--triodata_keep_pheno`
File with a list of participants to keep in analysis. If this file is provided together with the file to `--triodata_fam` argument, triodata_define process is skipped and these two files are used instead.
Example argument:
```
--triodata_keep_pheno s3://lifebit-featured-datasets/projects/gel/siteqc/keep.txt
```
Example file content:
```
sample_1 sample_1
sample_2 sample_2
sample_3 sample_3
HG002 HG002
HG003 HG003 
HG004 HG004
```
#### `--triodata_fam`
A .fam file describing families of participants selected to keep in analysis. If this file is provided together with the file to `--triodata_keep_pheno` argument, triodata_define process is skipped and these two files are used instead.
Example argument:
```
--triodata_fam s3://lifebit-featured-datasets/projects/gel/siteqc/example.fam
```
Example file content:
```
family_1 sample_1 sample_2 sample_3 1 -99
family_1 sample_2 0 0 1 -99
family_1 sample_3 0 0 2 -99
family_2 HG002 HG003 HG004 2 -99 
family_2 HG003 0 0 1 -99
family_2 HG004 0 0 2 -99
```
## Optional SiteQC pipeline arguments
### `--updfile`
File with a list of participants' platekeys to be exluded from the analysis.
Optional input to triodata_define process.
File should have same structure as file for `--included_samples`.
We were not provided with any example file forr that argument.

### `--mend_err_p3_keep_fam`
A .fam file of families that should be kept after mend_err_* filtering processes. By default such file is generated by mend_dist process, but user can provide a different file to be used instead. In this case mend_dist and mend_err_p2 processes will be skipped.
Example argument:
```
--triodata_fam s3://lifebit-featured-datasets/projects/gel/siteqc/MendelFamilies_4SD.fam
```
Example file content:
```
FID
family_1
family_2
```
### `--xx_sample_ids`
Path to file that contains the XX sample ids. Must be single column with sample ids, without a header.
Example argument:
```
--xx_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xx.txt
```
### `--xy_sample_ids`
Path to file that contains the XY sample ids. Must be single column with sample ids, without a header.
Example argument:
```
--xy_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xy.txt
```

### Path to 1000 genomes files.
SiteQC aggregate annotation step also includes some data from 1000Genomes files of corresponding chromosomes. Depending on the input files, only needed chromosomes will be fetched from the path. Path to s3 files can be different as well as their names, for this reason the argument should contain both full path and also name pattern of the files where only chromosome name would vary. For the current version of pipeline we are using files located on 1000Genomes bucket. The path is split into two parts - before and after the chromosome name in file name. Chromosome name is inserted programmatically by nextflow and whole paths are created during pipeline run.

Default path to 1000G files:
```
--s3_path_1kg_start 's3://1000genomes/release/20130502/ALL.'
--s3_path_1kg_end '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
```

## Other SiteQC parameters
The following parameters are not intended to be changed by users, but can be changed if needed.
### `--king`
Boolean value, set to `T` by default. This activates some part of final aggregate_annotation process R script. For current pipeline implementation should be always set to true.
### bcftools queries and awk commands
These arguments are virtually outsourced parts of scripts to make code cleaner. Changing any of those without deep knowledge will most probably break the pipeline.
| parameter | Default value |
| --------- | ------------- |
| `query_format_start` | `'%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n'` |
| `query_format_miss1` | `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n'` |
| `query_include_miss1` | `'GT="./." & FORMAT/DP=0'` |
| `query_format_miss2` |  `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n'` |
| `query_exclude_miss2` | `'GT~"\\."'` |
| `query_exclude_median_gq` | `'GT~"\\."'` |
| `query_exclude_med_cov_nonmiss` |  `'GT~"\\."'` |
| `query_format_med_cov_all` | `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n'` |
| `awk_expr_med_cov_all` | `'{sum=0; n=split(\$0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a); median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2; print \$1"\t"median}'` |
| `query_format_med_cov_nonmisss` | `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n'` |
| `awk_expr_med_cov_nonmiss` | `'{sum=0; n=split(\$0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a); median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2; print \$1"\t"median}'` |
| `query_format_median_gq` | `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n'` |
| `awk_expr_median_gq` | `'{sum=0; n=split(\$0,a)-1; for(i=2;i<=n;i++) sum+=a[i]; asort(a); median=n%2?a[n/2+1]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99}; print \$1"\t"median}'` |
| `query_format_ab_ratio_p1` | `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n'` |
| `query_include_ab_ratio_1` | `'GT="het" & binom(FMT/AD) > 0.01'` |
| `query_format_ab_ratio_p2`| `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n'` |
| `query_include_ab_ratio_2` |  `'GT="het"'` |
| `query_format_pull_ac` | `'%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC %INFO/AC \n'` |
| `awk_expr_miss1` | `'BEGIN{FS=" "} {print $1" "NF-1}'` |
| `awk_expr_miss2` | `'BEGIN{FS=" "} {print $1" "NF-1}'` |
| `awk_expr_ab_ratio_1` | `'{print $1"\t"NF -1}'` |
| `awk_expr_pull_1kg_p3_sites` | `'/^##/ {next} { print $1"\t"$2"\t"$4"\t"$5}'` |
| `mend_err_p1_rset_missing_var_ids` | `'@:#,\\\$r,\\\$a'` |
| `mend_err_p1_vcf_half_call` | `'m'` |
| `mend_err_p1_new_id_max_allele_len` | `'60 missing'` |


## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/siteqc`](https://hub.docker.com/r/nfcore/siteqc/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/siteqc`](https://hub.docker.com/r/nfcore/siteqc/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
