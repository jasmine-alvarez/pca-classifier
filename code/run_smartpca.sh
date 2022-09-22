#!/bin/bash

# pre-check
infile=data/acs_mini_project.vcf.gz ; bcftools view -m2 -M2 -v snps $infile | \
bcftools +prune -l 0.1 | bcftools +fill-tags | bcftools filter -e 'AF < 0.01' | \
bcftools filter -e 'F_MISSING >= 0.05' | bgzip -c > data/acs_mini_project.filtered.vcf.gz

plink --bfile data/acs_mini_project --pca  --snps-only --remove data/unlabeled_inds.txt--chr-set 22 no-xy

bcftools +fixref data/acs_mini_project.filtered.vcf.gz -- -d -f data/GRCh37.fa.gz -i All_20151104.vcf.g
plink --vcf data/acs_mini_project.filtered.vcf.gz --remove data/unlabeled_inds.txt --snps-only --recodeA --out data/acs_mini_project.2plink
plink --vcf data/acs_mini_project.filtered.vcf.gz --make-bed --snps-only --remove data/unlabeled_inds.txt --out data/acs_mini_project.2plink
# check

for i in {1..22} ; do \
shapeit -check -V data/acs_mini_project.chr${i}.unphased.vcf.gz \
-M data/1000GP_Phase3/genetic_map_chr${i}_combined_b37.txt \
-R data/1000GP_Phase3/1000GP_Phase3_chr${i}.hap.gz \
data/1000GP_Phase3/1000GP_Phase3_chr${i}.legend.gz \
data/1000GP_Phase3/1000GP_Phase3.sample \
--output-log data/acs_mini_project.chr${i}.alignments ; done

plink --vcf data/acs_mini_project.filtered.vcf.gz --make-bed --out data/acs_mini_project ; \
cp data/acs_mini_project.bim data/acs_mini_project.pedsnp ; \
cp data/acs_mini_project.fam data/acs_mini_project.pedind

plink --bfile data/acs_mini_project --pca --chr-set 2 no-xy

while read chr start stop fname; do \
echo ${fname} ; comm -12  <(cut -f2 -d' ' output/${fname}.imputed.gen | sort | grep -v \\.) \
<(cut -f2 -d' ' data/working/acs_mini_project.${fname}.phased.haps | sort | grep -v \\.) > ${fname}.to.filter ; \
plink --data output/${fname}.imputed --make-bed --snps-only --extract ${fname}.to.filter --oxford-single-chr ${chr} --out data/${fname} ; \
done < data/acs_mini_project.AIMS.txt

touch allfiles.txt ; \
while read chr start stop fname; do \
echo "data/${fname}.bed data/${fname}.bim data/${fname}.fam" >> allfiles.txt ; \
done < data/acs_mini_project.AIMS.txt

cp data/${fname}.bim data/{fname}.pedsnp ; \
cp data/${fname}.fam data/{fname}.pedind

fname=acs_mini_project ; \
plink --bfile data/${fname} --recodeA --out data/${fname}

# remember to remove NA from data
#
plink --vcf data/acs_mini_project.filtered.vcf.gz --make-bed --snps-only --remove data/unlabeled_inds.txt --out data/acs_mini_project.2plink

INFILE="data/acs_mini_project.2plink" ; \
../EIG/bin/smartpca.perl -i ${INFILE}.bed \
-a ${INFILE}.pedsnp -b ${INFILE}.pedind \
-p ${INFILE}.plot -o ${INFILE}.pca \
-e ${INFILE}.eval -l ${INFILE}.log
