#!/bin/bash


plink --file data/hapmap1 --make-bed --out data/hapmap1

../EIG/bin/smartpca.perl -i data/hapmap1.bed -a data/hapmap1.pedsnp -b data/hapmap1.pedind -o data/hapmap1.pca -p data/hapmap1.plot -e data/hapmap1.eval -l data/hapmap1.log
