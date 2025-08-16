# PanK-Pipeline: A Pan-genome _K_-mer Pipeline for Population Analysis

## Installation
The pipeline [PanK-Pipeline](https://github.com/caixu0518/PanK-Pipeline) is installation-free but requires dependencies: 

Required:
1. [jellyfish](https://github.com/gmarcais/Jellyfish). In the present pipeline,[jellyfish](https://github.com/gmarcais/Jellyfish) is mainly used to quickly generate _k_-mers from resequencing reads and perform _k_-mer query.

Optional：
1. [plink](https://www.cog-genomics.org/plink2/) (v1.90b6.21). [plink](https://www.cog-genomics.org/plink2/) is mainly used when performing population structure and PCA analysis with _k_-mers.
2. [VCF2Dis](https://doi.org/10.1093/gigascience/giaf032)(VCF2Dis-1.54). [VCF2Dis](https://doi.org/10.1093/gigascience/giaf032) is used to make the phylogenetic tree based on _k_-mer presence and absence matrix.
3. [faststructure](https://github.com/rajanil/fastStructure). [faststructure](https://github.com/rajanil/fastStructure) is used to make the population structure analysis. The present pipeline recomeneded a docker repository (dockerbiotools/faststructure).   

## Inputs

## Outputs


## Introduction
### Step1. Workflow for constructing polymorphic k-mers using 30 B. rapa genome assemblies.


### Step2: Workflow for identifying representative k-mers across the B. rapa species.


### Step3：Application of species-representative k-mers for population structure analysis in B. rapa.



## Citations


