

#!/bin/bash

#Wget Clytia genome fasta gz and gff3 gz (and protein sequences nr)
wget --content-disposition https://www.dropbox.com/sh/ek9j50haemzv57y/AAChG0nQjsfoP9jnWvwDDf56a?dl=0

unzip hirise

gunzip clytia_hic_output.fasta.gz

#Make indexed genome
/home/tchari/samtools-1.9/samtools faidx /home/tchari/clytiaRef/clytia_hic_output.fasta

gunzip transcripts.fa.transdecoder.genome.gff3.gz

python makeTransdecoderGTF.py #Make gtf and bed file for genes

gzip transcripts.fa.transdecoder.genome_geneAndTrans.gtf
gzip clytia_hic_output.fasta

kb ref -i new_trans_index.idx -g new_trans_t2g.txt -f1 new_trans_cdna.fa -f2 new_trans_intron.fa -c1 new_trans_cdna_t2c.txt -c2 new_trans_intron_t2c.txt --overwrite --workflow lamanno \
clytia_hic_output.fasta.gz \
transcripts.fa.transdecoder.genome_geneAndTrans.gtf.gz
#Marimba_merged_transcript_models.gff3.gz  

gunzip clytia_hic_output.fasta.gz

/home/tchari/bedtools getfasta -fi clytia_hic_output.fasta -bed transdecoder_genes.bed -name > new_trans_clytia_all_genes.fa

cat new_trans_clytia_all_genes.fa | cut -f1 -d":" > new_trans_clytia_all_genes_correct.fa

/home/tchari/bedtools subtract -a transdecoder_genes.bed -b transdecoder_exons.bed > transdecoder_introns.bed
/home/tchari/bedtools getfasta -fi clytia_hic_output.fasta -bed transdecoder_introns.bed -name > new_trans_clytia_all_introns.fa

cat new_trans_clytia_all_introns.fa | cut -f1 -d":" > new_trans_clytia_all_introns_correct.fa


/home/tchari/bedtools getfasta -fi clytia_hic_output.fasta -bed transdecoder_exons.bed -name > new_trans_clytia_all_exons.fa
cat new_trans_clytia_all_exons.fa | cut -f1 -d":" > new_trans_clytia_all_exons_correct.fa

