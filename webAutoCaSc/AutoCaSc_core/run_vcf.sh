#!/bin/sh

for input in /home/johann/Bernt_VCF/pedigrees/*
do
filename="${input##*/}"
pedigree_name="${filename##*_}"
family="${pedigree_name%.*}"
echo "$family"
python AutoCaSc_core/AutoCaSc_vcf.py score_vcf -v /home/johann/Bernt_VCF/bernt_vcf_coding -p '/home/johann/Bernt_VCF/pedigrees/pedigree_'"$family"'.ped' -g /home/johann/tools/slivar/gnotate/gnomad.hg38.v2.zip -j /home/johann/PycharmProjects/AutoCaSc_project_folder/AutoCaSc_core/data/slivar-functions.js -o '/home/johann/Bernt_VCF/scored/'"$family"'_scored.csv' -a GRCh38 -s /home/johann/tools/slivar/slivar
done


