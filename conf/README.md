# SiteQC Profiles and Tests

## General description
Detailed nf-core guide for nextflow configuartion files can be found [here](https://github.com/nf-core/configs).

### Config files
Normal Nextflow pipeline execution by default takes into account config file `nextflow.config` if it is present in the base dir of the pipeline.
Additional configs can be specified via command line argument `-c my_config.config`, or defined in main `nextflow.config` file as following:
```
includeConfig 'conf/aws.config'
```
All configs should be kept at `conf/` folder of the pipeline repo, which you are browsing now.
### Profiles
Profile is a mode in which pipeline is executed, and allows to run pipeline with different set of paramters/configs by changing just the profile name in the command line:
```
nextflow run main.nf -profile my_profile
```
Such profiles should be defined in `nextflow.config` in scope `profiles` ([docs](https://www.nextflow.io/docs/latest/config.html#config-profiles)).
```
profiles {
  standard { includeConfig 'conf/aws.config' }
  test { includeConfig 'conf/test.config' }
}
```
By default, CloudOS uses profile called `standard` to run nextflow pipelines on the platform, so for each pipeline to be run o CloudOS it is preferrable to have such a profile with a corresponding `aws.config` file. The latter is a conventional name for a config used to provide parameters for running Nextflow pipeline over AWS clod compute system from CloudOS.

Additionally, any custom profile can be added. Typical use case would be to add one or multiple test profiles that are using predefined demo data to test pipeline viability. If there are more than one way to run the pipeline, a test profile for each use case should be created. These test profiles are used to run github CI tests for submitted commits and PRs. The tests are defined in [ci.yml](https://github.com/lifebit-ai/siteqc/blob/master/.github/workflows/ci.yml) file of the repository. 


## SiteQC test profiles
SiteQC pipeline has several modes of execution based on the input files provided. Below you can find definition of each test profile.
All these profiles are included in CI testing that is automatically triggered on each commit or PR submission.
 
### Test 1 - `test_full`
**Profile name and aliases:** `test`, `test1`, `test_full`.

**Description:**
Main expected way to run the pipeline. User should provide all required inputs as defined in documentation.

The test can be run as following:
```
nextflow run main.nf -profile test
```
Which is equivalent to:
```
nextflow run main.nf --input s3://lifebit-featured-datasets/projects/gel/siteqc/input.csv \
                     --triofile s3://lifebit-featured-datasets/projects/gel/siteqc/extended_interpreted_trios_grch38.tsv \
                     --included_samples s3://lifebit-featured-datasets/projects/gel/siteqc/sampleList.txt \
                     --xx_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xx.txt \
                     --xy_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xy.txt
```
The first run will take ~4 minutes, use `-resume` option in subsequent runs to save time on not running unchanged parts of pipeline.

## Test 2 - `test_triodata_keep_provided`
**Profile name and aliases:** `test2`, `test_triodata_keep_provided`.

**Description:**

Run the pipeline with minimum processes possible by providing needed files for two conditional executions:
- provide triodata `keep_pheno` and `family` files to skip `triodata_define` process;
- provide `MendelFamilies_4SD.fam` file to skip `mend_err_p2` and `mend_dist` processes.
```
nextflow run main.nf -profile test2 -resume
```
Which is equivalent to:
```
nextflow run main.nf --input s3://lifebit-featured-datasets/projects/gel/siteqc/input.csv \
                     --included_samples s3://lifebit-featured-datasets/projects/gel/siteqc/sampleList.txt \
                     --xx_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xx.txt \
                     --xy_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xy.txt \
                     --triodata_keep_pheno s3://lifebit-featured-datasets/projects/gel/siteqc/keep.txt \
                     --triodata_fam s3://lifebit-featured-datasets/projects/gel/siteqc/example.fam \
                     --mend_err_p3_keep_fam s3://lifebit-featured-datasets/projects/gel/siteqc/MendelFamilies_4SD.fam \
                     -resume #optional, to save time by getting the cached big 1kg big file.
```
## Test 3 - `test_triodata_provided`

**Profile name and aliases:** `test3`, `test_triodata_provided`.

**Description:**

Run the pipeline as in test 1, but with skipping `triodata_define` process only, by providing triodata `keep_pheno` and `family` files.
```
nextflow run main.nf -profile test3 -resume
```
Which is equivalent to:
```
nextflow run main.nf --input s3://lifebit-featured-datasets/projects/gel/siteqc/input.csv \
                     --included_samples s3://lifebit-featured-datasets/projects/gel/siteqc/sampleList.txt \
                     --xx_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xx.txt \
                     --xy_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xy.txt \
                     --triodata_keep_pheno s3://lifebit-featured-datasets/projects/gel/siteqc/keep.txt \
                     --triodata_fam s3://lifebit-featured-datasets/projects/gel/siteqc/example.fam \
                     -resume #optional
```


## Test 4 - `test_keep_provided`

**Profile name and aliases:** `test4`, `test_keep_provided`.

**Description:**

Run the pipeline as in test 1, but with skipping `mend_err_p2` and `mend_dist` processes, by providing `MendelFamilies_4SD.fam` file.
```
nextflow run main.nf -profile test4 -resume
```
Which is equivalent to:
```
nextflow run main.nf --input s3://lifebit-featured-datasets/projects/gel/siteqc/input.csv \
                     --triofile s3://lifebit-featured-datasets/projects/gel/siteqc/extended_interpreted_trios_grch38.tsv \
                     --included_samples s3://lifebit-featured-datasets/projects/gel/siteqc/sampleList.txt \
                     --xx_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xx.txt \
                     --xy_sample_ids s3://lifebit-featured-datasets/projects/gel/siteqc/xy.txt \
                     --mend_err_p3_keep_fam s3://lifebit-featured-datasets/projects/gel/siteqc/MendelFamilies_4SD.fam \
                     -resume #optional
```
