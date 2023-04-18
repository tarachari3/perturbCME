#!/bin/bash
# generate count matrices

kb count --verbose \
-i /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/index.idx \
-g /home/ggorin/ref/refdata-gex-GRCh38-2020-A/t2g_grch38.txt \
-x 10xv3 \
-o ../counts/jiang_drugcombo/L1/ \
-t 30 -m 30G \
-c1 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/cdna_t2c.txt \
-c2 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L1_S13_L002_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L1_S13_L002_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L1_S13_L003_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L1_S13_L003_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L1_S13_L004_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L1_S13_L004_R2_001.fastq.gz & \

kb count --verbose \
-i /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/index.idx \
-g /home/ggorin/ref/refdata-gex-GRCh38-2020-A/t2g_grch38.txt \
-x 10xv3 \
-o ../counts/jiang_drugcombo/L2/ \
-t 30 -m 30G \
-c1 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/cdna_t2c.txt \
-c2 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L2_S14_L002_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L2_S14_L002_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L2_S14_L003_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L2_S14_L003_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L2_S14_L004_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L2_S14_L004_R2_001.fastq.gz & \

kb count --verbose \
-i /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/index.idx \
-g /home/ggorin/ref/refdata-gex-GRCh38-2020-A/t2g_grch38.txt \
-x 10xv3 \
-o ../counts/jiang_drugcombo/L3/ \
-t 30 -m 30G \
-c1 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/cdna_t2c.txt \
-c2 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L3_S15_L002_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L3_S15_L002_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L3_S15_L003_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L3_S15_L003_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L3_S15_L004_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L3_S15_L004_R2_001.fastq.gz & \

kb count --verbose \
-i /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/index.idx \
-g /home/ggorin/ref/refdata-gex-GRCh38-2020-A/t2g_grch38.txt \
-x 10xv3 \
-o ../counts/jiang_drugcombo/L4/ \
-t 30 -m 30G \
-c1 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/cdna_t2c.txt \
-c2 /home/ggorin/ref/refdata-gex-GRCh38-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L4_S16_L002_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L4_S16_L002_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L4_S16_L003_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L4_S16_L003_R2_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L4_S16_L004_R1_001.fastq.gz \
../counts/jiang_drugcombo/LocalFolder/SCMT01/MULT-19-R1-L4_S16_L004_R2_001.fastq.gz
