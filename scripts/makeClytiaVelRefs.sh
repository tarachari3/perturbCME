
#!/bin/bash

#Wget Clytia genome fasta gz
wget --content-disposition https://data.caltech.edu/tindfiles/serve/130f876e-f813-4cd6-bcfd-a41a8244cc18/


#Wget Clytia gtf transcriptome (add gtf to CaltechData) gz doi 20053
#wget --content-disposition https://data.caltech.edu/tindfiles/serve/4a7c1a67-e640-43f3-a2ef-c4b83dd49e69/
#Update with Marimba gtf doi: 20055
wget --content-disposition https://data.caltech.edu/tindfiles/serve/46362e66-626d-4b24-88db-aa1f1a7e01ee/

kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --overwrite --workflow lamanno \
genome.fa.gz \
Marimba_merged_transcript_models.gtf.gz  


