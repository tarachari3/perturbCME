
import scipy.io as sio
import pandas as pd
import numpy as np
from scipy.sparse import vstack

import scanpy as sc
import anndata
import loompy

samples = ['FT-SA168S1S4','FT-SA168S5S8']


data_path = "../counts/clytia_starv/"
out_path = "../counts/clytia_starv/loom/"

meta_path = "../metadata/"


#Read in S/U mtx and concatenate

#'/counts_unfiltered', 'spliced.mtx' ,'unspliced.mtx', 'spliced.genes.txt'
spliced =[]
unspliced = []

sNames = []
uNames = []

#Read in scanpy metadata
adata = anndata.read_h5ad(meta_path+'fedStarved_withUMAPPaga.h5ad') #https://data.caltech.edu/records/mm6y6-g4569/files/fedStarved_withUMAPPaga.h5ad.gz?download=1

#Filter for Fed only
adata = adata[adata.obs['fed']=='True']

barcodes = list(adata.obs_names) #Fed barcodes

for i in range(len(samples)):
	samp = samples[i]

	s = sio.mmread(data_path+samp+'/counts_unfiltered/spliced.mtx')

	u = sio.mmread(data_path+samp+'/counts_unfiltered/unspliced.mtx')

	origSNames = [x+'-'+str(i+1) for x in list(pd.read_csv(data_path+samp+'/counts_unfiltered/spliced.barcodes.txt',header=None)[0])]
	origUNames = [x+'-'+str(i+1) for x in list(pd.read_csv(data_path+samp+'/counts_unfiltered/unspliced.barcodes.txt',header=None)[0])]

	origSNames_sub = list(set(origSNames).intersection(origUNames))
	origUNames_sub = list(set(origSNames).intersection(origUNames))

	sToKeep = [i in barcodes for i in origSNames_sub]
	uToKeep = [i in barcodes for i in origUNames_sub]


	s = s.tocsr()[sToKeep,:]
	spliced += [s.todense()]

	u = u.tocsr()[uToKeep,:]
	unspliced += [u.todense()]


	sNames += np.array(origSNames_sub)[sToKeep].tolist()
	uNames += np.array(origUNames_sub)[uToKeep].tolist()

ds = loompy.connect(data_path+samp+'/counts_unfiltered/adata.loom')
geneNames = ds.ra['gene_id'] #np.array(list(pd.read_csv(data_path+samp+'/counts_unfiltered/spliced.genes.txt',header=None)[0]))
ds.close()

allS = vstack(spliced).toarray()
allU = vstack(unspliced).toarray()

print('Matrix sizes:')
print(allS.shape)
print(allU.shape)

print('Barcodes:')
print(len(sNames))
print(len(uNames))





#Save loom prior to gene filtering
fname = out_path+'clytia_Fed_allGenes.loom'


retAdata = anndata.AnnData(
	X=allS,
	layers={
		'spliced': allS,
		'unspliced': allU
	},
	obs=pd.DataFrame({'barcode': np.array(sNames)},index=np.array(sNames)),
	var=pd.DataFrame({'gene_name': geneNames},index=geneNames)
)

retAdata.write_loom(fname)
























