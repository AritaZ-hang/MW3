## STEP01 SNP CALLING

#!/bin/bash

reference=~/Mus_musculus.GRCm38.88.fasta
name=lung_wt_tumor_100wbin_pseudobulks

cd `pwd`

## sort and index
samtools sort -o $name"_sorted.bam" $name".bam"
samtools index $name"_sorted.bam"

# bcftools vcf calling
bcftools mpileup -Ou -f $reference $name"_sorted.bam" | bcftools call -vm -Oz > "./"$name".vcf.gz"

gunzip $name".vcf.gz"

## STEP 02  Divide the VCF file

#!/bin/bash

cd `pwd`

vcf_sep=~/vcf_sep.py
name=lung_wt_tumor_100wbin_pseudobulks
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

name=lung_wt_tumor_100wbin_pseudobulks
tissue_bam=~/$name"_sorted.bam"
bclist=~/$name"_bc.tsv"
vcf_path=~/vcf
vartrix=~/vartrix
reference=~/Mus_musculus.GRCm38.88.fasta

## processing on the regular chromosomes of mus musculus
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
cbn=~/Cbn_matrix.R
name=lung_wt_tumor_100wbin_pseudobulks

dir_path=~/
vcf_path=~/vcf
bclist=~/$name"_bc.tsv"

Rscript $cbn $dir_path $vcf_path $bclist $samplename
