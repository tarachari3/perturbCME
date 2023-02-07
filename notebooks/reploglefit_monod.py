import sys
sys.path.insert(0, '/home/ggorin/projects/monod/src/')

import monod
from monod import preprocess, extract_data, cme_toolbox, inference, analysis

import pandas as pd
import numpy as np
import loompy as lp
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import scipy.stats

meta_path = "/home/tchari/metadata/"
meta = pd.read_csv(meta_path+'replogle_GSM4367984_concat_identities.csv')


dataset_meta = ['allcrispr']
print('dataset_meta: ', dataset_meta)
print()

subcluster_names = [['sgPMVKa_sgFDPSb'],['sgPMVKa_sgNegCtrl3b'],['sgNegCtrl02093a_sgFDPSb'],
					['sgHUS1_2a_sgFDPS_2b'],['sgHUS1_2a_sgNegCtrl3b'],['sgNegCtrl02093a_sgFDPS_2b'],
					['sgNegCtrl02093a_sgNegCtrl3b'],
					['sgNegCtrl4a_sgNegCtrl3b'],['sgNegCtrl8a_sgNegCtrl3b'],
					['sgNegCtrl5a_sgNegCtrl3b'],['sgNegCtrl1a_sgNegCtrl3b'],
					['sgNegCtrl9a_sgNegCtrl3b'],['sgNegCtrl4a_sgNegCtrl3b','sgNegCtrl8a_sgNegCtrl3b',
					'sgNegCtrl5a_sgNegCtrl3b','sgNegCtrl1a_sgNegCtrl3b','sgNegCtrl9a_sgNegCtrl3b']]

sub_names_only = ['_'.join(n) for n in subcluster_names]

cluster_names = []
dataset_names = ['replogle_'+dataset_meta[0]+'_'+y  for y in sub_names_only]   #To save
print('dataset_names: ', dataset_names)
print('len(dataset_names): ',len(dataset_names))
print()

raw_data_locations = [dataset_meta[0] for y in sub_names_only]
transcriptome_filepath = '/home/tchari/perturbCME/notebooks/gg_200525_genome_polyA_cum_3'

spliced_layer = 'spliced'
unspliced_layer = 'unspliced'
gene_attr = 'gene_name'
cell_attr = 'barcode'

attribute_names=[(unspliced_layer,spliced_layer),gene_attr,cell_attr]

loom_filepaths = ['/home/tchari/counts/replogle_crispr/loom/'+x+'.loom' for x in raw_data_locations] 
print('loom_filepaths: ',loom_filepaths)

n_datasets = len(loom_filepaths)


cf = []
thr_lb = [1e4]*len(dataset_meta)*2

fig1,ax1 = plt.subplots(1,len(dataset_meta)*2,figsize=(15,5))

for k in range(len(dataset_meta)):
	filename = loom_filepaths[len(subcluster_names)*k]
	dataset_name = raw_data_locations[len(subcluster_names)*k]
	
	with lp.connect(filename,mode='r') as ds:
		S = ds.layers[spliced_layer][:]
		U = ds.layers[unspliced_layer][:]
		gene_names = ds.ra[gene_attr]
		bcs = ds.ca[cell_attr]
		n_cells = S.shape[1]
		monod.preprocess.knee_plot(S+U,ax1[k],viz=True,thr=thr_lb[k])
		cf_ = ((S+U).sum(0)>thr_lb[k])
		
		n_annot_bcs = meta['cell_barcode'].sum()
		annot_bcs_in_loom = meta['cell_barcode'].isin(bcs).sum()
		annot_bcs_in_filt_loom = meta['cell_barcode'].isin(bcs[cf_]).sum()
		print(f'Dataset {dataset_name}. \n\t{len(bcs)} barcodes in loom, {cf_.sum()} pass filter. {n_annot_bcs} in annotations; of these, {annot_bcs_in_loom} in loom and {annot_bcs_in_filt_loom} in filtered loom.')
		
		#if k==0:
		for subcluster in subcluster_names:
			annot_bcs = meta[(meta['guide_identity'].isin(subcluster))]['cell_barcode']
			cf.append(np.isin(bcs,annot_bcs) & cf_)
			print(f'\t{subcluster}: {len(annot_bcs)} cells in annotations. {np.isin(bcs,annot_bcs).sum()} in loom. {cf[-1].sum()} pass filter.')



	ax1[k].set_title(dataset_name)
	

fig_dir = './figs/'
fig_string = fig_dir + 'kneeplots_all_replogle.png'
fig1.tight_layout()
plt.savefig(fig_string,dpi=450)


print("Cells per subcluster: ", [x.sum() for x in cf])
print()


import logging, sys
logging.basicConfig(stream=sys.stdout)
log = logging.getLogger()
log.setLevel(logging.INFO)


dir_string,dataset_strings = monod.preprocess.construct_batch(loom_filepaths, transcriptome_filepath,  dataset_names, attribute_names=attribute_names, batch_location='./fits',meta='replogle_crispr', batch_id=1, n_genes=3000,exp_filter_threshold=None,cf=cf)

print(dir_string)
print('DATASET STRINGS')
print()
print(dataset_strings)

#Define bounds
phys_lb = [-2.0, -1.8, -1.8 ] #-1.0, -1.8, -1.8
phys_ub = [4.2, 2.5, 2.5] #4.2, 2.5, 3.5
samp_lb = [-9, -4] #-7.5, -2
samp_ub = [-4, 1.5] #-5.5, 0
# gridsize = [5,6]
gridsize = [20,21]

result_strings = []
for i in range(n_datasets):
	fitmodel = monod.cme_toolbox.CMEModel('Bursty','Poisson')
	inference_parameters = monod.inference.InferenceParameters(phys_lb,phys_ub,samp_lb,samp_ub,gridsize,\
				dataset_strings[i],fitmodel,use_lengths = True,
				gradient_params = {'max_iterations':20,'init_pattern':'moments','num_restarts':5})
	search_data = monod.extract_data.extract_data(loom_filepaths[i], transcriptome_filepath, dataset_names[i],
				dataset_strings[i], dir_string, dataset_attr_names=attribute_names,cf=cf[i])
	full_result_string = inference_parameters.fit_all_grid_points(60,search_data)

	result_strings.append(full_result_string)


print(result_strings)



#Currently showing results for old run
for i in range(n_datasets):
	sr = monod.analysis.load_search_results(result_strings[i])
	sd = monod.analysis.load_search_data(dir_string+'/'+dataset_names[i]+'/raw.sd')
	fig1,ax1 = plt.subplots(1,1)
	sr.find_sampling_optimum()
	sr.plot_landscape(ax1)

	fig1,ax1 = plt.subplots(1,1)
	sr.plot_KL(ax1)

	sr.plot_gene_distributions(sd,marg='joint')

	_=sr.chisquare_testing(sd,threshold=1e-3)
	sr.resample_opt_viz()
	sr.resample_opt_mc_viz()
	sr.chisq_best_param_correction(sd,Ntries=4,viz=False,threshold=1e-3) 

	sr.compute_sigma(sd,num_cores=60)
	sr.plot_param_L_dep(plot_errorbars=False,plot_fit=True)
	sr.plot_param_marg()
	
	monod.analysis.make_batch_analysis_dir([sr],dir_string)
	sr.update_on_disk()






