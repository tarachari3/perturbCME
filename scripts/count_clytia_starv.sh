#!/bin/bash
# generate count matrices

kb count --verbose \
-i /home/tchari/clytiaRef/index.idx \
-g /home/tchari/clytiaRef/t2g.txt \
-x 10xv2 \
-o ../counts/clytia_starv/FT-SA168S1S4/ \
-t 30 -m 30G \
-w /home/tchari/clytiaRef/10xv2_whitelist.txt \
-c1 /home/tchari/clytiaRef/cdna_t2c.txt \
-c2 /home/tchari/clytiaRef/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/clytia_starv/FT-SA16888_S1_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16888_S1_L004_R2_001.fastq.gz \
../counts/clytia_starv/FT-SA16889_S2_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16889_S2_L004_R2_001.fastq.gz \
../counts/clytia_starv/FT-SA16890_S3_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16890_S3_L004_R2_001.fastq.gz \
../counts/clytia_starv/FT-SA16891_S4_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16891_S4_L004_R2_001.fastq.gz & \

kb count --verbose \
-i /home/tchari/clytiaRef/index.idx \
-g /home/tchari/clytiaRef/t2g.txt \
-x 10xv2 \
-o ../counts/clytia_starv/FT-SA168S5S8/ \
-t 30 -m 30G \
-w /home/tchari/clytiaRef/10xv2_whitelist.txt \
-c1 /home/tchari/clytiaRef/cdna_t2c.txt \
-c2 /home/tchari/clytiaRef/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/clytia_starv/FT-SA16892_S5_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16892_S5_L004_R2_001.fastq.gz \
../counts/clytia_starv/FT-SA16893_S6_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16893_S6_L004_R2_001.fastq.gz \
../counts/clytia_starv/FT-SA16894_S7_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16894_S7_L004_R2_001.fastq.gz \
../counts/clytia_starv/FT-SA16895_S8_L004_R1_001.fastq.gz \
../counts/clytia_starv/FT-SA16895_S8_L004_R2_001.fastq.gz
