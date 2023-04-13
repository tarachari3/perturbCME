from datetime import date
import time
t1 = time.time()
import warnings
warnings.filterwarnings("ignore")

# Can add this function later to seq_cme_inference
# Write input txt file

def genRunInput(fname='example_input.txt',dataDir = './',outDir='./',loomName = '',tranName = '',polyaA='15',tInd='0',
				filt='0.01, 0.01, 350, 350, 3, 3',exclude='',attList="[['spliced','unspliced','Gene','Barcode']]",
				seed='42', nGenes='5000', ind='0', gList='', lb='-2, -1.8, -1.8', ub='4.2, 2.5, 2.5',
				restart='1',init='moments',length='True',niter='20',nCu='10',nlambda='11',
				lbSamp='-9, -4',ubSamp='-4, 1.5',suffix='1',creator='tc',ncor='20',override=''):
	'''
	Generate run file for CME inference procedure. 
	Parameters listed in order of appearance in run text file.
	For a dry-run set ind = '-1'.

	Parameters:
	fname: filepath for run file
	dataDir: filepath for data directory
	outDir: filepath for output directory
	loomName: loom filenames
	tranName: transcriptome filepath
	polyaA: min polyA count
	tInd: column of transcriptome file (lengths if 0, polyA count if 1)
	filt: min threshold for mean, max threshold for max, mean threshold for max; odd is U, even is S
	exclude:  result files (with genes to exclude)
	attList: list of lists for loom attribute names
	seed: gene selection seed
	nGenes: number of genes to select
	ind: loom file to analyze
	gList: result files to use for gene names
	lb: log10 of lower bound on burst size, splice rate, degradation rate
	ub: upper bound on on burst size, splice rate, degradation rate
	restart: how many times to run the search for each gene
	init: method to start search; 'random' if nothing
	length: use gene length for Poisson sampling
	niter: max iterations of gradient descent
	nCu: number of C_u points to evaluate
	nlambda: number of lambda_s points to evalute
	lbSamp: lower limits of C_u and lambda_s
	ubSamp: upper limits of C_u and lambda_s
	suffix: folder directory suffix
	creator: directory creator name
	ncor: number of cores to use
	override: if empty, use today's date from computer, else use string provided

	Returns:
	Run file (txt)
	'''

	with open(fname, 'w') as f:
		f.write('#Parameter input for CME inference ')
		f.write('\n')
		f.write('dataset_directory : '+dataDir+' : folder with dataset loom files')
		f.write('\n')
		f.write('result_directory : '+outDir+' : where to put the result folder')
		f.write('\n')
		f.write('loom_filenames : '+loomName+' : filenames to integrate')
		f.write('\n')
		f.write('transcriptome_filename : '+tranName+' : transcriptome location')
		f.write('\n')
		f.write('polyA_threshold : '+polyaA+' : minimum polyA count to use for sampling function')
		f.write('\n')
		f.write('transcriptome_ind : '+tInd+' : column of transcriptome file to use (lengths if 0, polyA count if 1)')
		f.write('\n')
		f.write('filter_param : '+filt+' : min threshold for mean, max threshold for max, mean threshold for max; odd is U, even is S')
		f.write('\n')
		f.write('all_prev_results : '+exclude+': result files with gene names to exclude')
		f.write('\n')
		f.write('attribute_names : '+attList+' : list or list of lists with loom attribute names')
		f.write('\n')
		f.write('gene_sel_seed : '+seed+' : gene selection seed')
		f.write('\n')
		f.write('n_gen : '+nGenes+' : number of genes to select')
		f.write('\n')
		f.write('IND : '+ind+' : loom_filename to analyze')
		f.write('\n')
		f.write('gene_list : '+gList+' : set of result files to import to define gene names')
		f.write('\n')
		f.write('phys_lb : '+lb+' : log10 of lower bound on burst size, splice rate, degradation rate')
		f.write('\n')
		f.write('phys_ub : '+ub+' : upper bound on same')
		f.write('\n')
		f.write('search_restarts : '+restart+' : how many times to run the search for each gene')
		f.write('\n')
		f.write("init_pattern : "+init+" : whether to start the search using method of moments estimate or not; 'random' if not")
		f.write('\n')
		f.write('use_lengths : '+length+' : whether the Poisson sampling for unspliced mRNA should depend on gene length')
		f.write('\n')
		f.write('maxiter : '+niter+' : number of iterations of gradient descent to perform')
		f.write('\n')
		f.write('n_pt1 : '+nCu+' : number of C_u points to evaluate')
		f.write('\n')
		f.write('n_pt2 : '+nlambda+' : number of lambda_s points to evalute')
		f.write('\n')
		f.write('samp_lb : '+lbSamp+' : lower limits of C_u and lambda_s')
		f.write('\n')
		f.write('samp_ub : '+ubSamp+' : upper limits of C_u and lambda_s')
		f.write('\n')
		f.write('ID_suffix : '+suffix+' : folder directory suffix')
		f.write('\n')
		f.write('creator : '+creator+' : directory creator name, can also be used for generic metadata')
		f.write('\n')
		f.write('NCOR : '+ncor+' : number of cores to use')
		f.write('\n')
		f.write("date_override : "+override+": if empty, use today's date from computer. if not, use the given string")




# a e s t h e t i c s

#all
col_gold = [203/255,197/255,149/255]
col_gray = [116/255,112/255,113/255]
col_red = [212/255,107/255,75/255]

### SI figures
ms = 0.7
lw = ms
ms_fail = 1
alf = 0.3
alf_fail = 0.5
sifig_aesth = (ms,lw,ms_fail,alf,alf_fail)

### Body figures
ms = 2
lw = ms
ms_fail = 2
alf = 0.3
alf_fail = 0.8
msfig_aesth = (ms,lw,ms_fail,alf,alf_fail)

markerstyle = None



def import_precomputed_(filestring):
	#Import result.pickle file
	with open(filestring,'rb') as f:
		precomp = pickle.load(f)
	return precomp

def check_at_bounds(result_data,phys_params,thr=0.01):
	#Check for genes where parameter fits are close to boundary of grid
	x = np.any(np.logical_or(
		result_data.search_params.lb_log+thr > phys_params,
		result_data.search_params.ub_log-thr < phys_params) ,1)
	return x

def get_name(fname,iden="_20x21"):
	#Get all result files for run
	start = fname.find(data_path+proj_fold+out_fold+creator+'_'+override+'_') + len(data_path+proj_fold+out_fold+creator+'_'+override+'_')
	end = fname.find(iden)
	substring = fname[start:end]
	
	return substring




from seq_cme_inference import *
from driver import *
from itertools import combinations
from matplotlib.colors import ListedColormap
import glob
import numpy as np
import pandas as pd
import seaborn as sns



data_path = '/home/tchari/counts/'

in_fold = 'loom/'
out_fold = 'loom_res/'

transcriptome = 'gg_200525_genome_polyA_cum_3'

#Dataset-dependent

# Only use neg control samples with >500 cells
control_looms =  ['crisprNegCtrl10_NegCtrl0__NegCtrl10_NegCtrl0','crisprNegCtrl11_NegCtrl0__NegCtrl11_NegCtrl0',
				 'crisprNegCtrl1_NegCtrl0__NegCtrl1_NegCtrl0','crisprNegCtrl0_NegCtrl0__NegCtrl0_NegCtrl0'] 

control_combo_loom = ['crisprNegCtrl10_NegCtrl0__NegCtrl10_NegCtrl0_NegCtrl11_NegCtrl0__NegCtrl11_NegCtrl0_NegCtrl1_NegCtrl0__NegCtrl1_NegCtrl0']
#STUDY uses sgNegCtrl4a_sgNegCtrl3b  only as control

perturb_looms = ['crisprDUSP9_ETS2__DUSP9_ETS2',
		 'crisprNegCtrl0_ETS2__NegCtrl0_ETS2_ETS2_NegCtrl0__ETS2_NegCtrl0','crisprDUSP9_NegCtrl0__DUSP9_NegCtrl0',
				'crisprCBL_NegCtrl0__CBL_NegCtrl0', 'crisprCBL_CNN1__CBL_CNN1',
				'crisprNegCtrl0_CNN1__NegCtrl0_CNN1_CNN1_NegCtrl0__CNN1_NegCtrl0']  #Maybe change to sgHUS1_2a_sgFDPS_2b, sgNegCtrl02093a_sgFDPS_2b

proj_fold = 'norman_crispr/'   

override='220517' #'220324'


creator = 'tc'




for i in range(0,len(control_looms)):
	genRunInput(fname=control_looms[i]+'_input.txt',dataDir =data_path+proj_fold+in_fold,
				outDir=data_path+proj_fold+out_fold,loomName = control_looms[i],
				nCu='20',nlambda='21', restart='5', niter='40', ncor='40',
				attList= "[['spliced','unspliced','Gene','Barcode']]", tranName = transcriptome,override=override) #gList= selGenes+','+filtGenes,
	

for i in range(0,len(control_looms)):
	inference_workflow(control_looms[i]+'_input.txt')



#nCu='1',nlambda='1',
#lbSamp='-7.333333333333333, -1.7999999999999998',ubSamp='-7.333333333333333, -1.7999999999999998',
for i in range(0,len(perturb_looms)):
	genRunInput(fname=perturb_looms[i]+'_input.txt',dataDir =data_path+proj_fold+in_fold,
				outDir=data_path+proj_fold+out_fold,loomName = perturb_looms[i], 
				nCu='20',nlambda='21', restart='5', niter='40', ncor='40',
				attList= "[['spliced','unspliced','Gene','Barcode']]", tranName = transcriptome,override=override) #gList= selGenes+','+filtGenes,

for i in range(0,len(perturb_looms)):
	inference_workflow(perturb_looms[i]+'_input.txt')
