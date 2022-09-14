#!/bin/bash

# TODO LD-pruned to r2 < 0.1, high callrate, biallelic, AF > 0.01
plink --file data/amp.prune --make-bed --out data/acs_mini_project 

# 1. filter to biallelic, common, high call rate
bcftools view -m2 -M2 -v snps data/acs_mini_project.vcf.gz | \
bcftools filter -e 'AF < 0.1' | \
bcftools filter -e 'F_MISSING > 0.05' | \
bgzip -c > data/acs_mini_project.fltr.vcf.gz

# 2. prune
plink --vcf data/acs_mini_project.fltr.vcf.gz  --indep 50 5 2
mv plink* data/

# 2a. convert back to vcf
plink --vcf data/acs_mini_project.fltr.vcf.gz --extract data/plink.prune.in \
--make-bed --keep-allele-order --out data/acs_mini_project.prun

plink --bfile data/acs_mini_project.prun --recode vcf --keep-allele-order --out data/acs_mini_project.tmp

grep -v ^23 data/acs_mini_project.tmp.vcf | \
bgzip -c > data/acs_mini_project.tmp.vcf.gz

# 3. normalize
#bcftools norm -f data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz | \
vt normalize data/acs_mini_project.tmp.vcf.gz -r data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz -n -o data/acs_mini_project.norm.vcf

plink --file data --recodeAD

../EIG/bin/smartpca.perl -i data/acs_mini_project.prun.bed -a data/acs_mini_project.prun.pedsnp -b data/acs_mini_project.prun.pedind -k 60 -t 60 -o data/acs_mini_project.pca -p data/acs_mini_project.plot -e data/acs_mini_project.eval -l data/acs_mini_project.log

printf %.10f\n 1 
