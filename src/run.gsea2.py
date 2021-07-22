import os, sys
from optparse import OptionParser
import argparse
import shutil

import gp
import kwat
from gp import data
from gsea import run_gsea

import pandas as pd

def main():
	usage="%prog [options]" + "\n"
	ap = argparse.ArgumentParser()
	ap.add_argument("-g","--gct",action="store",dest="GCT",help="GCT file.")
	ap.add_argument("-c","--cls",action="store",dest="CLS",help="CLS file.")
	ap.add_argument("-m","--gmt",action="store",dest="GMT",help="GMT file.")
	options = ap.parse_args()

	os.mkdir("gsea_results")

	sc_ge_sa = gp.data.GCT(options.GCT) # Parse GCT file
	cls_file = gp.data.CLS(options.CLS) # Parse CLS file
	gm = kwat.gmt.read([options.GMT]) # Parse GMT file
	
	sc_ge_sa.index=sc_ge_sa.index.droplevel(1) # Drop gene descriptions
	sa_la = pd.Series(cls_file.class_assignments) # extract class assignments from CLS object
	sa_la.index = sc_ge_sa.columns # Assign sample names to classes
	
	nu_se_st = run_gsea(
        sc_ge_sa,  # Gene-by-sample score; DataFrame
        sa_la,  # Sample label; Series
        gm,  # Gene sets; set-to-genes dict or DataFrame
        no="-0-",  # Normalization method; "-0-", "0-1", "1234", "log"
        ra="ic",  # Ranking method; "ic", "si", "co", "tt", "di", "ra", "lo"
        mi=5,  # Minimum gene set size; int
        ma=500,  # Maximum gene set size; int
        we=1.0,  # Weight used for "ks" and "auc"; float
        me="js",  # Enrichment method; "ks", "auc", "js"
        se=1729,  # Random seed; int
        pe="label",  # Permutation type; "gene_set", "label"
        n_pe=1000,  # Number of permutations; int
        n_pl=25,  # Number of extreme gene sets to plot; int
        ad=None,  # Additional gene sets to plot; list of str
        pa="gsea_results", # directory path to write the gene-set-by-statistic and plots; str
	)
	
	nu_se_st

if __name__ == '__main__':
	main()
