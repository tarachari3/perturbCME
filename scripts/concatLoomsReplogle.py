
import scipy.io as sio
import pandas as pd
import numpy as np
from scipy.sparse import vstack

import scanpy as sc
import anndata
import loompy

samples = ['lane1','lane2','lane3','lane4','lane5','lane6']


data_path = "../counts/replogle_crispr/"
out_path = "../counts/replogle_crispr/loom/"

meta_path = "../metadata/"


#Read in S/U mtx and concatenate

#'/counts_unfiltered', 'spliced.mtx' ,'unspliced.mtx', 'spliced.genes.txt'
spliced =[]
unspliced = []

sNames = []
uNames = []

for i in range(len(samples)):
	samp = samples[i]

	s = sio.mmread(data_path+samp+'/counts_unfiltered/spliced.mtx')
	spliced += [s]

	u = sio.mmread(data_path+samp+'/counts_unfiltered/unspliced.mtx')
	unspliced += [u]

	sNames += [x+'-'+str(i+1) for x in list(pd.read_csv(data_path+samp+'/counts_unfiltered/spliced.barcodes.txt',header=None)[0])]
	uNames += [x+'-'+str(i+1) for x in list(pd.read_csv(data_path+samp+'/counts_unfiltered/unspliced.barcodes.txt',header=None)[0])]

#geneNames = np.array(list(pd.read_csv(data_path+samp+'/counts_unfiltered/spliced.genes.txt',header=None)[0]))
ds = loompy.connect(data_path+samp+'/counts_unfiltered/adata.loom')
geneNames = ds.ra['gene_name'] #np.array(list(pd.read_csv(data_path+samp+'/counts_unfiltered/spliced.genes.txt',header=None)[0]))
ds.close()

allS = vstack(spliced) # .toarray()
allU = vstack(unspliced) #.toarray()

print('Matrix sizes:')
print(allS.shape)
print(allU.shape)

print('Barcodes:')
print(len(sNames))
print(len(uNames))


#Read in metadata

meta = pd.read_csv(meta_path+'GSM4367984_exp6.cell_identities.csv')
meta['guide_identity'] = [i+'|'+j for i,j in zip(meta.A, meta.B)] #Make global identity strings/labels

#pair = ['NegCtrl10_NegCtrl0__NegCtrl10_NegCtrl0','NegCtrl11_NegCtrl0__NegCtrl11_NegCtrl0','NegCtrl1_NegCtrl0__NegCtrl1_NegCtrl0']
pair = []

remain = np.unique(meta.guide_identity)
remain = [[i] for i in remain if i not in pair]

#All conditions/paired conditions
assigns = remain #+ [pair]

#Comment out to not split up controls
#assigns = [['NegCtrl10_NegCtrl0__NegCtrl10_NegCtrl0'],['NegCtrl11_NegCtrl0__NegCtrl11_NegCtrl0'],['NegCtrl1_NegCtrl0__NegCtrl1_NegCtrl0']]

#For each drug condition get cell barcodes/counts and save loom file
for a in assigns:

	barcodes = list(meta['cell_barcode'][meta['guide_identity'].isin(a)])

	sfilt = [sNames.index(x) for x in barcodes]
	print('Condition: ', a)
	print('Sp. Filt:', len(sfilt))
	ufilt = [uNames.index(x) for x in barcodes]


	subS = allS.tocsr()[sfilt,:]
	subU = allU.tocsr()[ufilt,:]

	print(subS.shape)
	print(subU.shape)

	if subS.shape[0] > 0:
		#Save loom files in data_path
		names = '_'.join(a)
		fname = out_path+'crispr'+names+'.loom'

		#row_attrs = { "Gene": geneNames } #genes
		#col_attrs = { "Barcode": np.array(barcodes) } #cells

		retAdata = anndata.AnnData(
			X=subS,
			layers={
				'spliced': subS,
				'unspliced': subU
			},
			obs=pd.DataFrame({'Barcode': np.array(barcodes)},index=np.array(barcodes)),
			var=pd.DataFrame({'Gene': np.array(geneNames)},index=np.array(geneNames))
		)

		retAdata.write_loom(fname)























