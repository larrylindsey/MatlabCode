#!/bin/bash

d1=/groups/chklovskii/home/nuneziglesiasj/Projects/em_denoising/data/10x10x10_cropped/segmentation_maps
d2=/groups/chklovskii/home/nuneziglesiasj/Projects/datatemp/tao/segment_maps
p1=1_115_lda11

for t in $@
do
    outfile=/groups/chklovskii/home/nuneziglesiasj/Projects/em_denoising/data/10x10x10_cropped/evals/bwr_$t.csv
    gt_fn=/groups/chklovskii/home/nuneziglesiasj/Projects/em_denoising/data/10x10x10_cropped/vtk_format/segmented.vtk
    ws_fns="{ '$d1/${p1}_12500.tif.ws',  '$d1/${p1}_25000.tif.ws', '$d1/${p1}_50000.tif.ws', '$d1/${p1}_100000.tif.ws', '$d1/${p1}_mito_18750.tif.ws', '$d1/${p1}_mito_37500.tif.ws', '$d1/${p1}_mito_75000.tif.ws', '$d1/${p1}_mito_150000.tif.ws', '$d2/mf3.ws', '$d2/mf3_no_mito.ws' }"
    amb_fns="{ '$d1/${p1}_12500.tif.amb_%03i.raw',  '$d1/${p1}_25000.tif.amb_%03i.raw', '$d1/${p1}_50000.tif.amb_%03i.raw', '$d1/${p1}_100000.tif.amb_%03i.raw', '$d1/${p1}_mito_18750.tif.amb_%03i.raw', '$d1/${p1}_mito_37500.tif.amb_%03i.raw', '$d1/${p1}_mito_75000.tif.amb_%03i.raw', '$d1/${p1}_mito_150000.tif.amb_%03i.raw', '$d2/mf3.ws.amb_%03i.raw', '$d2/mf3_no_mito.ws.amb_%03i.raw' }"
    ts=$t
    index_sets="{ {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115}, {1:500, 1:500, 1:115} }"
    colname_base="{ 'lda12k', 'lda25k', 'lda50k', 'lda100k', 'ldam18k', 'ldam37k', 'ldam75k', 'ldam150k', 'mf3', 'mf3m' }"
    trim=10
    min_body_size=3000
    markers_per_body=20
    submatrix_multiplier=5
    d=3


    #echo qsub -N bwr$t -cwd -V -b y -pe batch 4 /groups/chklovskii/home/nuneziglesiasj/Projects/em_recon/code/mcc/all_bwr_cluster \""$outfile"\" \""$gt_fn"\" \""$ws_fns"\" \""$amb_fns"\" \""$ts"\" \""$index_sets"\" \""$colname_base"\" \""$trim"\" \""$min_body_size"\" \""$markers_per_body"\" \""$submatrix_multiplier"\" \""$d"\"
    qsub -N bwr$t -cwd -V -b y -pe batch 1 /groups/chklovskii/home/nuneziglesiasj/Projects/em_recon/code/mcc/all_bwr_cluster \""$outfile"\" \""$gt_fn"\" \""$ws_fns"\" \""$amb_fns"\" \""$ts"\" \""$index_sets"\" \""$colname_base"\" \""$trim"\" \""$min_body_size"\" \""$markers_per_body"\" \""$submatrix_multiplier"\" \""$d"\"
done
