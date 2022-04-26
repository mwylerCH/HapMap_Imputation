
# Introduction

Various SNP imputation are available. Most methods require a VCF input format, knowledge about the pedigree or the linkage map. The advantage of *fastPHASE* (Scheet & Stephens, 2006), is the purely genotyping and position based method.  
*fastPHASE* is a Hidden Markov model based imputation method that clusters the haplotypes around the missing genotyping and imputes based on similarity. 

# Quick start

Download and prepare *HapMap_Imputation*:
```
git clone https://github.com/mwylerCH/HapMap_Imputation/
chmod +x HapMap_Imputation/HapMap_Imputation.pl
export PATH="$PATH:$HOME/HapMap_Imputation"

```
Usage 
```
HapMap_Imputation.pl FILE.hapmap > FILE.imputed
```


# Method

## software details

The new script requires a genotyping information in a hapmap format (comma separated, http://augustogarcia.me/statgen-esalq/Hapmap-and-VCF-formats-and-its-integration-with-onemap/#hapmap). 

In a first step, *HapMap_Imputation* counts the occourence of each nucleotide at every single genotyped position. The most common nucleotide is defined as major allele, the second is defined as minor allele. Missing genotyping information is excluded. In the case major and minor alleles accour at the same number, the nucleotide of the reference *cs10* (available as GCF_900626175.1 on NCBI) is choosen as major allele.

Subsequently, *HapMap_Imputation* sorts markers by position and parses the hapmap into the required *fastPHASE* input format (requirements explained in the *fastPHASE* manual). Briefly, *HapMap_Imputation* splits the haplotypes into two separate rows, converts major and minor allesles into 0 and 1 respectively and produces temporary files for each chromosomes.

During the third step, *HapMap_Imputation* downloads the latest *fastPHASE* version and runs the imputation using 8 cores in parallel. *fastPHASE* is runned with ten random starts of the imputation algorithm.

After imputation, *HapMap_Imputation* reverse the 0 and 1 coding into the major and minor nucleotide, respectively. Subsequently, the two haplotypes are combined and the the separate chromosome merged.



## Software testing

To test the capabilities of *HapMap_Imputation*, we used hapmap files from two populations ( a total of 431 individuals genotyped with 5536 markers). We subsequently juxtaposed the imputed and the raw files by comparing each single position for each individual.

# Results

The used files contained 26.42% missing values (630'438 of 2'386'016 of marker X individual combinations). Out of the 630'438 missing genotypes, *HapMap_Imputation* was able to impute 565'415 SNPs (89.7% of the missing values). The remaining 65'023 (10% of the missing values) were not imputable. This is mostly because of the limited genotyping over all individuals of the specific SNPs. A small fraction 0.27% (6'466) of all SNPs contains an incompatibilitiy with the original raw data. 


# Reference
Scheet, Paul, and Matthew Stephens. "A fast and flexible statistical model for large-scale population genotype data: applications to inferring missing genotypes and haplotypic phase." The American Journal of Human Genetics 78.4 (2006): 629-644.
