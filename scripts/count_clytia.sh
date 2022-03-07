#!/bin/bash
# generate count matrices

kb count --verbose \
-i /home/tchari/clytiaRef/index.idx \
-g /home/tchari/clytiaRef/t2g.txt \
-x 10xv3 \
-o ../counts/clytia_stim/FT-SA22418/ \
-t 30 -m 30G \
-w /home/tchari/clytiaRef/10xv3_whitelist.txt \
-c1 /home/tchari/clytiaRef/cdna_t2c.txt \
-c2 /home/tchari/clytiaRef/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/clytia_stim/FT-SA22418_S1_L001_R1_001.fastq.gz \
../counts/clytia_stim/FT-SA22418_S1_L001_R2_001.fastq.gz \
../counts/clytia_stim/FT-SA22418_S1_L002_R1_001.fastq.gz \
../counts/clytia_stim/FT-SA22418_S1_L002_R2_001.fastq.gz  & \


kb count --verbose \
-i /home/tchari/clytiaRef/index.idx \
-g /home/tchari/clytiaRef/t2g.txt \
-x 10xv3 \
-o ../counts/clytia_stim/FT-SA22419/ \
-t 30 -m 30G \
-w /home/tchari/clytiaRef/10xv3_whitelist.txt \
-c1 /home/tchari/clytiaRef/cdna_t2c.txt \
-c2 /home/tchari/clytiaRef/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../counts/clytia_stim/FT-SA22419_S2_L001_R1_001.fastq.gz \
../counts/clytia_stim/FT-SA22419_S2_L001_R2_001.fastq.gz \
../counts/clytia_stim/FT-SA22419_S2_L002_R1_001.fastq.gz \
../counts/clytia_stim/FT-SA22419_S2_L002_R2_001.fastq.gz

