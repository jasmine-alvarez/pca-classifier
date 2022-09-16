#!/bin/bash

# 0. phasing
for i in {1..22} ; do bcftools view -r ${i} data/acs_mini_project.vcf.gz | \
bcftools view -m2 -M2 -v snps | bgzip -c > data/acs_mini_project.chr${i}.vcf.gz ; done

for i in {1..22} ; do plink --vcf data/acs_mini_project.chr${i}.vcf.gz \
--missing --out data/chr${i}.missing ; done

for i in {1..22} ; do echo $i ; \
bcftools +fill-tags data/acs_mini_project.chr${i}.vcf.gz | \
bcftools view -e 'F_MISSING >= 0.05' -S data/acs_mini_project.chr${i}.samples | \
bgzip -c > data/acs_mini_project.chr${i}.unphased.vcf.gz ; done

for i in {1..22} ; do \
shapeit -check -V data/acs_mini_project.chr${i}.unphased.vcf.gz \
-M data/1000GP_Phase3/genetic_map_chr${i}_combined_b37.txt \
-R data/1000GP_Phase3/1000GP_Phase3_chr${i}.hap.gz \
data/1000GP_Phase3/1000GP_Phase3_chr${i}.legend.gz \
data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log acs_mini_project.chr${i}.alignments ; done

for i in {1..22} ; do \
shapeit -V data/acs_mini_project.chr${i}.unphased.vcf.gz \
-M data/1000GP_Phase3/genetic_map_chr${i}_combined_b37.txt \
-R data/1000GP_Phase3/1000GP_Phase3_chr${i}.hap.gz \
data/1000GP_Phase3/1000GP_Phase3_chr${i}.legend.gz \
data/1000GP_Phase3/1000GP_Phase3.sample \
--thread 8 -O data/acs_mini_project.chr${i}.phased \
--exclude-snp data/acs_mini_project.chr${i}.alignments.snp.strand.exclude ; done

#  1. filter
bcftools view -m2 -M2 -v snps data/acs_mini_project.vcf.gz | \
bcftools filter -e 'AF < 0.1' | \
bcftools filter -e 'F_MISSING > 0.05' | \
bgzip -c > data/acs_mini_project.fltr.vcf.gz

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
