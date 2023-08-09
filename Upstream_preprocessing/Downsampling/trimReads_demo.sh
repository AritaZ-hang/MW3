#!/bin/bash

scripts=~/grep_bam_usingReads.py
workdir=~

cd $workdir

grep_in_bam=${workdir}/mouse.bam
grep_out_bam=${workdir}/mouse_grep.bam

python $scripts $grep_in_bam $grep_out_bam

split_bam=~/splitbam.py
trim_in_bam=$grep_out_bam

samtools view -b -q 10 $trim_in_bam | samtools sort -> star_gene_exon_tagged_mapq10.bam

mkdir trimreads_15000
cd trimreads_15000

sample_name=$(basename `pwd`)
dropseq_root=~/Drop-seq_tools-2.5.1/
tmpdir=tmp
tmpdir2=tmp2

BAM_FILE=${workdir}/star_gene_exon_tagged_mapq10.bam
trimbc_root=~
str1="mouse_trim"
str2="reads"
file=.txt
file2=.bam
file3=_sorted.bam
file4=_sample.sam
file5=_sample.bam
file6=_dge.txt.gz
file7=_out_cell_readcounts.txt
file8=_dge.summary.txt

i=15000
mkdir tmp
samtools view -@ 16 -H $BAM_FILE > SAM_header
samtools view -@ 16 $BAM_FILE|LC_ALL=C grep -F -f $trimbc_root$str1$i$str2$file > filtered_SAM_body
cat SAM_header filtered_SAM_body > filtered.sam
samtools view -@ 16 -b filtered.sam > $str1$i$str2$file2 && rm filtered.sam filtered_SAM_body
samtools sort -t XC $str1$i$str2$file2 -o $str1$i$str2$file3 && rm $str1$i$str2$file2
python $split_bam $str1$i$str2$file3 ${workdir}/trimreads_15000/tmp/

for k in 1000 5000 10000 15000
do
    cd $tmpdir
    mkdir tmp2
    for j in *.bam
    do
        samtools view -@ 16 $j | shuf -n $k > $tmpdir2/filtered_SAM_body
        cat ../SAM_header $tmpdir2/filtered_SAM_body > $tmpdir2/$j$file4
        samtools view -@ 16 -b $tmpdir2/$j$file4 > $tmpdir2/$j$file5 && rm $tmpdir2/filtered_SAM_body $tmpdir2/$j$file4
        samtools sort -@ 16 $tmpdir2/$j$file5 > $tmpdir2/$j$file3 && rm $tmpdir2/$j$file5
        samtools index $tmpdir2/$j$file3
    done

    samtools merge ../$str1$k$file5 $tmpdir2/*_sorted.bam && rm -r $tmpdir2

    cd ../

    ${dropseq_root}/DigitalExpression -m 16g \
    I=$str1$k$file5 \
    CELL_BARCODE_TAG=XC \
    MOLECULAR_BARCODE_TAG=XM \
    O=$str1$k$file6 \
    SUMMARY=$str1$k$file8 \
    NUM_CORE_BARCODES=8000 \
    LOCUS_FUNCTION_LIST=INTRONIC \
    TMP_DIR=.

    ${dropseq_root}/BamTagHistogram I=$str1$k$file5 O=$str1$k$file7 TAG=XC
done