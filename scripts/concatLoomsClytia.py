
import scipy.io as sio
import pandas as pd
import numpy as np
from scipy.sparse import vstack

import scanpy as sc
import anndata
import loompy

samples = ['FT-SA22418','FT-SA22419']


data_path = "../counts/clytia_stim/"
out_path = "../counts/clytia_stim/loom/"

meta_path = "../metadata/"


#Read in S/U mtx and concatenate

#'/counts_unfiltered', 'spliced.mtx' ,'unspliced.mtx', 'spliced.genes.txt'
spliced =[]
unspliced = []

sNames = []
uNames = []

#Read in scanpy metadata
adata = anndata.read_h5ad(meta_path+'bus_stim.h5ad') #D1.1821

#Filter for SW only
adata = adata[adata.obs['condition'].isin(['SW'])] #Control cells only

barcodes = list(adata.obs_names) #SW barcodes

for i in range(len(samples)):
	samp = samples[i]

	s = sio.mmread(data_path+samp+'/counts_unfiltered/spliced.mtx')

	u = sio.mmread(data_path+samp+'/counts_unfiltered/unspliced.mtx')

	origSNames = [x+'-'+str(i+1) for x in list(pd.read_csv(data_path+samp+'/counts_unfiltered/spliced.barcodes.txt',header=None)[0])]
	origUNames = [x+'-'+str(i+1) for x in list(pd.read_csv(data_path+samp+'/counts_unfiltered/unspliced.barcodes.txt',header=None)[0])]

	sToKeep = [i in barcodes for i in origSNames]
	uToKeep = [i in barcodes for i in origUNames]

	s = s.tocsr()[sToKeep,:]
	spliced += [s.todense()]

	u = u.tocsr()[uToKeep,:]
	unspliced += [u.todense()]


	sNames += np.array(origSNames)[sToKeep].tolist()
	uNames += np.array(origUNames)[uToKeep].tolist()

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




#Make loom with SW no top 2000 hvgs (to fit sampling parameters)


#Get raw counts
adata_raw =  anndata.read_h5ad(meta_path+'jelly4stim_bus_combo_raw.h5ad') #D1.1814
sc.pp.filter_cells(adata_raw, min_counts=1)
sc.pp.filter_genes(adata_raw, min_counts=1)
#Get top hvgs
sc.pp.normalize_per_cell(adata_raw, counts_per_cell_after=1e4) 
sc.pp.log1p(adata_raw)

#
filter_result = sc.pp.filter_genes_dispersion(adata_raw.X, min_mean=0.0125, max_mean=4.5, min_disp=0.2)
#sc.pp.highly_variable_genes(adata_raw,n_top_genes=2000)
print(adata_raw)


#genes_to_keep = adata_raw[:,adata_raw.var['highly_variable']].var_names #actually want to remove these
genes_to_keep = adata_raw[:,filter_result.gene_subset].var_names #actually want to remove these
print(genes_to_keep[0])
print("XLOC_012650 in hvgs: ", 'XLOC_012650' in genes_to_keep)


sfilt = [sNames.index(x) for x in barcodes]
print(len(barcodes))
print('Sp. Filt:', len(sfilt))
ufilt = [uNames.index(x) for x in barcodes]


subS = allS[sfilt,:]
subU = allU[ufilt,:]

#Save loom prior to gene filtering
fname = out_path+'clytia_SWall_allGenes.loom'

#row_attrs = { "Gene": geneNames } #genes
#col_attrs = { "Barcode": np.array(barcodes) } #cells

retAdata = anndata.AnnData(
	X=subS,
	layers={
		'spliced': subS,
		'unspliced': subU
	},
	obs=pd.DataFrame({'Barcode': np.array(barcodes)},index=np.array(barcodes)),
	var=pd.DataFrame({'Gene': geneNames},index=geneNames)
)

retAdata.write_loom(fname)

#Filter out HVGs

to_keep = [i not in genes_to_keep for i in geneNames] # Genes not in hvgs
geneNamesSub = np.array(geneNames)[to_keep]

subS = subS[:,to_keep]
subU = subU[:,to_keep]

print('S Shape after filt: ', subS.shape)
print('U Shape after filt: ', subU.shape)

if subS.shape[0] > 0:
	#Save loom files in data_path

	fname = out_path+'clytia_SWall.loom'

	#row_attrs = { "Gene": geneNames } #genes
	#col_attrs = { "Barcode": np.array(barcodes) } #cells

	retAdata = anndata.AnnData(
		X=subS,
		layers={
			'spliced': subS,
			'unspliced': subU
		},
		obs=pd.DataFrame({'Barcode': np.array(barcodes)},index=np.array(barcodes)),
		var=pd.DataFrame({'Gene': geneNamesSub},index=geneNamesSub)
	)

	retAdata.write_loom(fname)


#Make looms with 'annos' cell types (broad categories: Neural, Nematocyte, etc) only (across all genes)


#Comment out if not splitting up controls
assigns = [[i] for i in np.unique(adata.obs['annos'])] #fix backslashes in name

#For each cell type get cell barcodes/counts and save loom file
for a in assigns:

	barcodes = list(adata[adata.obs['annos'].isin(a)].obs_names)

	sfilt = [sNames.index(x) for x in barcodes]
	print('Condition: ', a)
	print('Sp. Filt:', len(sfilt))
	ufilt = [uNames.index(x) for x in barcodes]


	subS = allS[sfilt,:]
	subU = allU[ufilt,:]

	print(subS.shape)
	print(subU.shape)

	if subS.shape[0] > 0:
		#Save loom files in data_path
		names = '_'.join(a)
		names = names.replace(' ','_')
		fname = out_path+'clytia'+names.replace('/','_')+'.loom'

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























