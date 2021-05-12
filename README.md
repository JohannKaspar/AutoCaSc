# AutoCaSc: evaluate variants in cases of NDD
If you use AutoCaSc, please cite the paper.
AutoCaSc is a command line tool for evaluating deleteriousness of genomic variants found in Trio Exome Sequencing in cases of neurodevelopmental disorders (NDD). Variants are assigned up to 15 points based on various parameters such as in silico predictions and inheritance patterns. The higher the score, the higher the probability that the variant is harmful. For detailed information on how the scores are determined, please refer to our paper.Variants are assigned up to 15 points based on various parameters such as in silico predictions and inheritance patterns. The higher the score, the higher the probability that the variant is harmful. For detailed information on how the scores are determined, please refer to our paper.

## Installation
1. Clone this repository
2. Install Python 3.8
3. Install pipenv
4. Navigate to the repository and run `pipenv install`

## Usage
There are different options for different usecases.
### Scoring single variants
If you have a set of variants that you want to quickly check or if you want to integrate AutoCaSc into a pipeline, use the command `single`. E.g. `pipenv run python AutoCaSc_core/AutoCaSc.py single --variant 1:55516888:G:GA --inheritance de_novo`
### Scoring multiple variants
If you want to score a list of variants, you can use the command `batch`. The command takes a table of variants in the form of comma-separated columns containing:
| Variant | Inheritance | Corresponding comphet variant | Family history |
| ----------- | ----------- | ----------- | ----------- |
| 14:93681589:G:GTT | comphet | 14:93673673:G:T | True |
| NM_015317.2:c.2216del | de_novo | | False |
### Scoring VCFs
`vcfAutoCaSc.py` can be used to score whole VCF files. These can be prefiltered using any tool of your choice. In this case inheritance patterns should be inserted into the INFO field of the variants. An interface to prefilter using [slivar](https://github.com/brentp/slivar) and [bedtools](https://bedtools.readthedocs.io/en/latest/index.html) is implemented in the script. Slivar needs to be installed seperately.
