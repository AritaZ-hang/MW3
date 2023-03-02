## STEP01 SNP CALLING

#!/bin/bash

reference=~/Mus_musculus.GRCm38.88.fasta
name=lung_wt_tumor
cd `pwd`

# bcftools vcf calling
bcftools mpileup -Ou -f $reference lung_wt_tumor.bam | bcftools call -vm -Oz > "./"$name".vcf.gz"

## STEP 02  Divide the VCF file

#!/bin/bash

cd `pwd`

vcf_sep=~/vcf_sep.py
tissue=lung_wt_tumor
outdir=~"/vcf"

if [ ! -d $outdir]; then
    mkdir $outdir
fi

python $vcf_sep -vp $tissue".vcf" -od $outdir -in False -out False


## STEP 03 Generating SNP by cell matrices for both reference allele and alternative alleles per chromosome

#!/bin/bash

cd `pwd`

outdir=chr
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

tissue=lung_wt_tumor
tissue_bam=~/$tissue".bam"
bclist=~/lung_bc.tsv
vcf_path=~/vcf
vartrix=~/vartrix
reference=~/Mus_musculus.GRCm38.88.fasta

## 
for chr in {1..19};
do
    mkdir chr/chr"$chr"_matrix
    $vartrix -v $vcf_path"/chr"$chr".vcf" --bam $tissue_bam -f $reference -c $bclist -s "coverage" --bam-tag CB
    mv out_matrix.mtx chr/chr"$chr"_matrix/out_matrix.mtx
    mv ref_matrix.mtx chr/chr"$chr"_matrix/ref_matrix.mtx
done

## STEP 04 Merge the matrices together

#!/bin/bash

cd `pwd`
cbn=/public/home/guogjgroup/ggj/rawdata/ZS/scripts/CNV/Cbn_matrix.R
samplename=lung_wt_tumor

dir_path=~/
vcf_path=~/vcf
bclist=~/lung_bc.tsv

Rscript $cbn $dir_path $vcf_path $bclist $samplename
