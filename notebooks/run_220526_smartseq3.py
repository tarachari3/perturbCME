from datetime import date
import time
import scipy.io as sio
import anndata
import loompy
import pandas as pd
import numpy as np
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

transcriptome = 'gg_200524_mouse_genome_polyA_cum_1'

#Dataset-dependent

# Only use neg control samples with >500 cells
smartseq_looms = ['smartSeq3Fibro']  #Maybe change to sgHUS1_2a_sgFDPS_2b, sgNegCtrl02093a_sgFDPS_2b

proj_fold = 'smartseq3_fibro/'   

override='220526' #'220324'


creator = 'tc'


genRunInput(fname=smartseq_looms[0]+'_input.txt',dataDir =data_path+proj_fold+in_fold,
			outDir=data_path+proj_fold+out_fold,loomName = smartseq_looms[0],nCu='20',nlambda='21',
			restart='5', niter='40',
			attList= "[['spliced','unspliced','gene_name','barcode']]", tranName = transcriptome,override=override)



inference_workflow(smartseq_looms[0]+'_input.txt')













