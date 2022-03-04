
#!/bin/bash

#Wget Clytia genome fasta gz
wget --content-disposition https://data.caltech.edu/tindfiles/serve/130f876e-f813-4cd6-bcfd-a41a8244cc18/


#Wget Clytia gtf transcriptome (add gtf to CaltechData) gz doi 20053
wget --content-disposition https://data.caltech.edu/tindfiles/serve/4a7c1a67-e640-43f3-a2ef-c4b83dd49e69/


kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno \
genome.fa.gz \
genes.gtf.gz 


