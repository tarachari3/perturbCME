#Make gtf for clytia with gene information for technical bias model

import pandas as pd
import numpy as np
import csv

#/home/tchari/gffread/gffread Marimba_merged_transcript_models.gff3 -T -o Marimba_merged_transcript_models.gtf
#python makeMarimbaGTF.py

gtf = pd.read_csv('Marimba_merged_transcript_models.gtf',sep='\t',header=None)

gff3 = pd.read_csv('Marimba_merged_transcript_models.gff3',sep='\t',header=None)

genes = gff3[gff3[2].isin(['gene'])] #Manually make gene entries gtf compatible

genes[8] = [i.replace('ID=','gene_id "') for i in genes[8]]
genes[8] = [i.replace(';Name=','"; gene_name "')+'"' for i in genes[8]]

concat = pd.concat([genes, gtf])

concat.to_csv('Marimba_merged_transcript_models_geneAndTrans.gtf',sep='\t',index=False,header=False,quoting=csv.QUOTE_NONE)

#Make bed file of genes only
gtf = pd.read_csv('Marimba_merged_transcript_models_geneAndTrans.gtf',header=None,sep='\t')
gtf = gtf[gtf[2].isin(['gene'])]
gtf[9] = [i[9:20] for i in gtf[8]] 
gtf[10] = gtf[[9,3,4]].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
gtf_sub = gtf[[0,3,4,10]]

gtf_sub.to_csv('Marimba_genes.bed',sep='\t',index=False,header=False)

