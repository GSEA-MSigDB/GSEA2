import os, sys
from optparse import OptionParser
import argparse
import shutil
import pickle

import gp
import kwat
from gp import data
from gsea import run_gsea

import pandas as pd
import numpy as np

def main():
	usage="%prog [options]" + "\n"
	ap = argparse.ArgumentParser()
	ap.add_argument("-g","--gct",action="store",dest="GCT",help="GCT file.")
	ap.add_argument("-c","--cls",action="store",dest="CLS",help="CLS file.")
	ap.add_argument("-m","--gmt",action="store",dest="GMT",help="GMT file.")
	ap.add_argument("-r","--rme",action="store",dest="rme",help="Ranking Metric.")
	ap.add_argument("-e","--eme",action="store",dest="eme",help="Enrichment Method.")
	ap.add_argument("-p","--perm",action="store",dest="perm",default="label",help="Permutation mode.")
	ap.add_argument("-n","--nperm",action="store",dest="nperm",default=1000,type=int,help="Number of permutations.")
	ap.add_argument("-a","--max",action="store",dest="max",default=500,type=int,help="Max gene set size.")
	ap.add_argument("-i","--min",action="store",dest="min",default=5,type=int,help="Min gene set size.")
	ap.add_argument("-w","--wgt",action="store",dest="weight",default=1.0,type=float,help="Weight for ks or auc enrichment method.")
	ap.add_argument("-l","--npl",action="store",dest="nplot",default=25,type=int,help="Number of enrichment results to plot.")
	ap.add_argument("-z","--rnd",action="store",dest="rnd",default=1729,type=int,help="Random seed used for permutations.")
	ap.add_argument("-j","--cpu",action="store",dest="cpu",default=1,type=int,help="Job CPU Count.")
	options = ap.parse_args()

	os.mkdir("gsea_results")

	sc_el_sa = gp.data.GCT(options.GCT) # Parse GCT file
	cls_file = gp.data.CLS(options.CLS) # Parse CLS file
	se_el_ = kwat.gmt.read([options.GMT]) # Parse GMT file
	
	sc_el_sa.index=sc_el_sa.index.droplevel(1) # Drop gene descriptions
	ta = pd.Series(cls_file.class_assignments) # extract class assignments from CLS object
	ta = np.asarray(ta) # convert to numpy array
	#ta.index = sc_el_sa.columns # Assign sample names to classes ## No longer works with numpy array class
	
	nu_se_st = run_gsea(
        ta,  # Sample label; Series
        sc_el_sa,  # Gene-by-sample score; DataFrame
        se_el_,  # Gene sets; set-to-genes dict or DataFrame
        fu=options.rme,  # Ranking method; "ic", "si", "co", "tt", "di", "ra", "lo"
        mi=options.min,  # Minimum gene set size; int
        ma=options.max,  # Maximum gene set size; int
        n_jo=options.cpu,
        we=options.weight,  # Weight used for "ks" and "auc"; float
        al=options.eme,  # Enrichment method; "ks", "auc", "js"
        pe=options.perm,  # Permutation type; "gene_set", "label"
        ra=options.rnd,  # Random seed; int
        n_pe=options.nperm,  # Number of permutations; int
        n_pl=options.nplot,  # Number of extreme gene sets to plot; int
        ad=None,  # Additional gene sets to plot; list of str
        pa="gsea_results", # directory path to write the gene-set-by-statistic and plots; str
	)

	pickle.dump(nu_se_st,open("gsea_results/complete_result.pkl", "wb"))

if __name__ == '__main__':
	main()


