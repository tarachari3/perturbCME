#!/bin/bash
# generate count matrices

kb count --verbose \
-i /home/ggorin/ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g /home/ggorin/ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x SMARTSEQ3 \
-o ../counts/smartseq3_fibro/kbOut/ \
-t 30 -m 30G \
-c1 /home/ggorin/ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 /home/ggorin/ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --overwrite --loom \
/home/tchari/counts/smartseq3_fibro/fastq/Smartseq3.Fibroblasts.NovaSeq.R1.fastq.gz  \
/home/tchari/counts/smartseq3_fibro/fastq/Smartseq3.Fibroblasts.NovaSeq.R2.fastq.gz

