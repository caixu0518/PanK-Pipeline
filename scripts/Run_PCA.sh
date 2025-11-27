#!/bin/bash

vcfIn=$1
prefix=$2

plink  --vcf  ${vcfIn}  --recode    --allow-extra-chr  --out  ${prefix}   --vcf-half-call  haploid

plink --file ${prefix}   --make-bed  --out ${prefix}  --allow-extra-chr

plink  --noweb --bfile  ${prefix}  --pca 20 --allow-extra-chr  --out plink.pca

#bgzip ${prefix}.vcf 

#tabix -p vcf ${prefix}.vcf.gz

#rm -rf ${prefix}.bed
#rm -rf ${prefix}.bim
#rm -rf ${prefix}.fam
#rm -rf ${prefix}.log
#rm -rf ${prefix}.map
#rm -rf ${prefix}.nosex
#rm -rf ${prefix}.ped

