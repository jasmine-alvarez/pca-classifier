#!/bin/bash

# TODO LD-pruned to r2 < 0.1, high callrate, biallelic, AF > 0.01
plink --file data/amp.prune --make-bed --out data/acs_mini_project 

bcftools view -m2 -M2 -v snps data/amp.prune.vcf.gz | \
bcftools filter -e 'AF < 0.1 & F_MISSING > 0.05' | \
bcftools norm -f data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz | \
bgzip -c > data/amp.prune.filtered.vcf.gz

../EIG/bin/smartpca.perl -i data/amp.prune.bed -a data/acs_mini_project.pedsnp -b data/acs_mini_project.pedind -k 60 -t 60 -o data/acs_mini_project.pca -p data/acs_mini_project.plot -e data/acs_mini_project.eval -l data/acs_mini_project.log

printf %.10f\n 1 
