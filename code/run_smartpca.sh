#!/bin/bash

# 1. filter to biallelic, common, high call rate
bcftools view -m2 -M2 -v snps data/acs_mini_project.vcf.gz | 
bcftools +fill-tags | bcftools filter -e 'AF <= 0.01' | \
bcftools filter -e 'F_MISSING > 0.05' | \
bgzip -c > data/acs_mini_project.no_miss.vcf.gz

# need to choose between bcftools vs plink pruning

plink --vcf data/acs_mini_project.no_miss.vcf.gz --make-bed --keep-allele-order --out data/acs_mini_project

cp data/acs_mini_project.bim data/acs_mini_project.pedsnp ; \
cp data/acs_mini_project.fam data/acs_mini_project.pedind

# remember to remove NA from data
../EIG/bin/smartpca.perl -i data/acs_mini_project.bed \
-a data/acs_mini_project.pedsnp -b data/acs_mini_project.pedind \
-k 200 -t 200 -p data/acs_mini_project.plot -o data/acs_mini_project.pca \
-e data/acs_mini_project.eval -l data/acs_mini_project.log

printf %.10f\n 1 
bcftools +setGT data/acs_mini_project.vcf.gz -- -t . -n 0
