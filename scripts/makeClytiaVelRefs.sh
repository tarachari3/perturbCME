
#!/bin/bash

#Wget Clytia genome fasta gz
wget --content-disposition https://data.caltech.edu/tindfiles/serve/130f876e-f813-4cd6-bcfd-a41a8244cc18/


#Wget Clytia gff3 transcriptome doi 1830
wget --content-disposition https://data.caltech.edu/tindfiles/serve/7299b2f2-a067-4ad6-9a0c-e92c73104841/

#Download genome.fa.fai

/home/tchari/gffread/gffread Marimba_merged_transcript_models.gff3 -T -o Marimba_merged_transcript_models.gtf

python makeMarimbaGTF.py #Make gtf and bed file for genes

gzip Marimba_merged_transcript_models_geneAndTrans.gtf 

kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --overwrite --workflow lamanno \
genome.fa.gz \
Marimba_merged_transcript_models_geneAndTrans.gtf.gz
#Marimba_merged_transcript_models.gff3.gz  

gunzip genome.fa.gz

/home/tchari/bedtools getfasta -fi genome.fa -bed Marimba_genes.bed -name > clytia_all_genes.fa

cat clytia_all_genes.fa | cut -f1 -d":" > clytia_all_genes_correct.fa

