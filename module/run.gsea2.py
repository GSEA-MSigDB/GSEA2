import os, sys, subprocess
from optparse import OptionParser
import argparse
import shutil
import json

import pandas
import numpy

def main():
	usage="%prog [options]" + "\n"
	ap = argparse.ArgumentParser()
	ap.add_argument("--libdir", action="store",
					dest="libdir", help="Working directory to load support library from.")
	ap.add_argument("--dataset",action="store",dest="dataset",help="Input Expression Dataset.")
	ap.add_argument("--gsdb",action="store",dest="gsdb",help="Gene Set Database File.")
	ap.add_argument("--nperm",action="store",dest="nperm",default=1000,type=int,help="Number of permutations.")
	ap.add_argument("--cls",action="store",dest="cls",help="CLS file.")
	ap.add_argument("--reverse",action="store",dest="reverse",default="False",help="Reverse the phenotype comparison defined in the CLS file.")
	ap.add_argument("--permute",action="store",dest="perm",default="False",help="Reverse the phenotype comparison defined in the CLS file.")
	ap.add_argument("--perm",action="store",dest="perm",default="sample",help="Permutation mode. Options are 'sample' (phenotype) and 'set' (gene set).")
	ap.add_argument("--collapse",action="store",dest="collapse",default="none",help="Method for computing mathematical collapse. Supports 'none' ('no collapse'), 'sum', 'mean', 'median', 'max', 'absmax'")
	ap.add_argument("--chip",action="store",dest="chip",default="none",help="Chip file used for performing collapse.")
	ap.add_argument("--metric",action="store",dest="rank_metric",help="Metric for ranking genes.")
	ap.add_argument("--method",action="store",dest="method",help="Enrichment Method. 'ks' (old GSEA) and 'js' (next gen GSEA) supported.")
	ap.add_argument("--weight",action="store",dest="weight",default=1.0,type=float,help="Weight for ks or auc enrichment method.")
	ap.add_argument("--max",action="store",dest="max",default=500,type=int,help="Max gene set size.")
	ap.add_argument("--min",action="store",dest="min",default=15,type=int,help="Min gene set size.")
	ap.add_argument("--seed",action="store",dest="seed",default=1729,type=int,help="Random seed used for permutations.")
	ap.add_argument("--nplot",action="store",dest="nplot",default=25,type=int,help="Number of enrichment results to plot.")
	ap.add_argument("--cpu",action="store",dest="cpu",default="1",type=int,help="Job CPU Count.")
	options = ap.parse_args()

	sys.path.insert(1, options.libdir)
	import GSEAlib

	os.mkdir("gsea_results")
	os.mkdir("gsea_results/input_files")


	## Parse GMT/GMX files from a list of inputs and create a name:members dict written out as a json file
	if options.gene_sets_db_list_filename != None:
		with open(options.gene_sets_db_list_filename) as f:
			gene_sets_dbfile_list = f.read().splitlines()

	genesets=GSEAlib.read_sets(gene_sets_dbfile_list)
	with open('gsea_results/input_files/set_to_genes.json', 'w') as path:
		json.dump(genesets, path,  indent=2)


	## Parse GCT file
	if options.input_gct_filename.split(".")[-1] == "gct":
		if options.collapse != "none":
			chip_file=read_chip(options.chip)
			input_ds = GSEAlib.collapse_dataset(options.input_gct_filename, chip_file, method=options.collapse)
		else:
			input_ds = GSEAlib.read_gct(options.input_gct_filename)
		input_ds=input_ds['data']
	else:
		input_ds=pandas.read_csv(options.input_gct_filename, sep='\t', index_col=0, skip_blank_lines=True)
		if "description" in input_ds.columns.str.lower():
			description_loc=input_ds.columns.str.lower().to_list().index('description')
			input_ds.drop(input_ds.columns[[description_loc]], axis = 1, inplace = True)
			input_ds.index.name="Name"
		if options.collapse != "none":
			chip_file=read_chip(options.chip)
			input_ds = GSEAlib.collapse_dataset(input_ds, chip_file, method=options.collapse)
			input_ds=input_ds['data']

	## Parse CLS file
	phenotypes=GSEAlib.read_cls(options.cls)
	phenotypes=GSEAlib.match_phenotypes(input_ds,phenotypes)
	if(options.reverse=="True"):
		phenotypes["Phenotypes"]=numpy.where((phenotypes["Phenotypes"]==0)|(phenotypes["Phenotypes"]==1), phenotypes["Phenotypes"]^1, phenotypes["Phenotypes"])
	phenotypes=phenotypes.sort_values('Phenotypes')

	## Order the dataset using the phenotypes and write out both files
	input_ds=input_ds.reindex(columns=phenotypes.index)
	input_ds.to_csv('gsea_results/input_files/gene_by_sample.tsv', sep = "\t")
	pandas.DataFrame(phenotypes['Phenotypes']).transpose().to_csv('gsea_results/input_files/target_by_sample.tsv',sep="\t", index=False)

	## Construct GSEA Settings json file
	gsea_settings={
		"number_of_permutations": options.nperm,
		"permutation": options.perm,
		"metric": options.rank_metric,
		"algorithm": options.method,
		"weight": options.weight,
		"maximum_gene_set_size": options.max,
		"minimum_gene_set_size": options.min,
		"random_seed": options.seed,
		"number_of_jobs": int(options.cpu),
		"number_of_extreme_gene_sets_to_plot": options.nplot,
		"gene_sets_to_plot": []
	}

	with open('gsea_results/input_files/gsea_settings.json', 'w') as path:
		json.dump(gsea_settings, path,  indent=2)

	## Run GSEA
	subprocess.check_output(['gsea', 'standard', 'gsea_results/input_files/gsea_settings.json', 'gsea_results/input_files/set_to_genes.json', 'gsea_results/input_files/target_by_sample.tsv', 'gsea_results/input_files/gene_by_sample.tsv', 'gsea_results'])

if __name__ == '__main__':
	main()
