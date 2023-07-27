import os
import sys
import subprocess
from optparse import OptionParser
from datetime import datetime
from zipfile import ZipFile
from os.path import basename
import argparse
import shutil
import json
import random
import pandas
import numpy
import dominate
from dominate.tags import *
from dominate.util import raw
from scipy.integrate import simps


# Better boolean command line parsing
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main():
    usage = "%prog [options]" + "\n"
    ap = argparse.ArgumentParser()
    ap.add_argument("--libdir", action="store",
                    dest="libdir", help="Working directory to load support library from.")
    ap.add_argument("--dataset", action="store", dest="dataset",
                    help="Input Expression Dataset.")
    ap.add_argument("--gsdb", action="store", dest="gsdb",
                    help="Gene Set Database File.")
    # ap.add_argument("--nperm", action="store", dest="nperm",
    #                 default=1000, type=int, help="Number of permutations.")
    # ap.add_argument("--cls", action="store", dest="cls", help="CLS file.")
    # ap.add_argument("--reverse", action="store", type=str2bool, nargs='?', const=True, dest="reverse",
    #                 default=False, help="Reverse the phenotype comparison defined in the CLS file.")
    # ap.add_argument("--perm", action="store", dest="perm", default="sample",
    #                 help="Permutation mode. Options are 'sample' (phenotype) and 'set' (gene set).")
    ap.add_argument("--collapse", action="store", dest="collapse", default="none",
                    help="Method for computing mathematical collapse. Supports 'none' ('no collapse'), 'sum', 'mean', 'median', 'max', 'absmax'")
    ap.add_argument("--chip", action="store", dest="chip",
                    default="none", help="Chip file used for performing collapse.")
    # ap.add_argument("--metric", action="store",
    #                 dest="rank_metric", help="Metric for ranking genes.")
    ap.add_argument("--alg", action="store", dest="method",
                    help="Enrichment Method. 'ks' (Kolmogorov-Smirnov, classic GSEA), 'ksa' (Kolmogorov-Smirnov area), kli, kliop, kliom")
    ap.add_argument("--exponent", action="store", dest="exponent", default=1.0,
                    type=float, help="Weight for ks or ksa enrichment method.")
    ap.add_argument("--max", action="store", dest="max",
                                    default=500, type=int, help="Max gene set size.")
    ap.add_argument("--min", action="store", dest="min",
                                    default=15, type=int, help="Min gene set size.")
    # ap.add_argument("--seed", action="store", dest="seed",
    #                 default="timestamp", help="Random seed used for permutations.")
    ap.add_argument("--ogllv", action="store", type=str2bool, nargs='?', const=True, dest="override",
                    default=False, help="Override reasonableness check for input dataset gene list size.")
    ap.add_argument("--nplot", action="store", dest="nplot", default=25,
                    type=int, help="Number of enrichment results to plot.")
    ap.add_argument("--zip", action="store", type=str2bool, nargs='?', const=True,
                                    dest="zip", default=False, help="Create ZIP bundle of results.")
    ap.add_argument("--cpu", action="store", dest="cpu",
                                    default=1, type=int, help="Job CPU Count.")
    options = ap.parse_args()

    sys.path.insert(1, options.libdir)
    import GSEAlib

    # Make a directory to store processed input files
    os.mkdir("input")

    # # Generate and set the random seed at the Python level and save it to pass to GSEA
    # if options.seed == "timestamp":
    #     options.seed = int(round(datetime.now().timestamp()))
    #     random.seed(options.seed)
    # else:
    #     options.seed = int(round(float(options.seed)))
    #     random.seed(options.seed)

    # Parse GCT file
    if options.dataset.split(".")[-1] == "gct":
        if options.collapse != "none":
            chip_file = GSEAlib.read_chip(options.chip)
            input_ds = GSEAlib.collapse_dataset(
                options.dataset, chip_file, method=options.collapse, drop=True)
            input_ds['mappings'].to_csv(
                'input/collapse_dataset_mapping_details.tsv', sep="\t", na_rep="No Symbol Mapping")
            collapse_length = input_ds['collapse_length']
        else:
            input_ds = GSEAlib.read_gct(options.dataset)
        input_length = input_ds['input_length']
        input_ds = input_ds['data']
    elif options.dataset.split(".")[-1] == "rnk":
        input_ds = pandas.read_csv(
            options.dataset, sep='\t', index_col=0, skip_blank_lines=True, header=None)
        if any(input_ds.index.str.startswith('#')):
            input_ds = input_ds.rename(columns=input_ds.iloc[int(
                numpy.where(input_ds.index.str.startswith('#'))[0])])
            input_ds = input_ds[input_ds.index.str.startswith('#') != True]
        else:
            input_ds = input_ds.rename(columns={1: "Preranked Metric"})
        input_ds.index.name = "Name"
        input_length = len(input_ds.index)
        if options.collapse != "none":
            chip_file = GSEAlib.read_chip(options.chip)
            input_ds = GSEAlib.collapse_dataset(
                input_ds, chip_file, method=options.collapse, drop=True)
            input_ds['mappings'].to_csv(
                'input/collapse_dataset_mapping_details.tsv', sep="\t", na_rep="No Symbol Mapping")
            input_length = input_ds['input_length']
            collapse_length = input_ds['collapse_length']
            input_ds = input_ds['data']
    else:
        input_ds = pandas.read_csv(
            options.dataset, sep='\t', index_col=0, skip_blank_lines=True)
        if "description" in input_ds.columns.str.lower():
            description_loc = input_ds.columns.str.lower().to_list().index('description')
            input_ds.drop(
                input_ds.columns[[description_loc]], axis=1, inplace=True)
        input_ds.index.name = "Name"
        input_length = len(input_ds.index)
        if options.collapse != "none":
            chip_file = GSEAlib.read_chip(options.chip)
            input_ds = GSEAlib.collapse_dataset(
                input_ds, chip_file, method=options.collapse, drop=True)
            input_ds['mappings'].to_csv(
                'input/collapse_dataset_mapping_details.tsv', sep="\t", na_rep="No Symbol Mapping")
            input_length = input_ds['input_length']
            collapse_length = input_ds['collapse_length']
            input_ds = input_ds['data']

    if len(input_ds) < 10000 and options.override == False and options.collapse == "none":
        sys.exit(print("Only ", len(input_ds), "genes were identified in the dataset.\nEither the dataset did not contain all expressed genes, or collapse dataset may need to be run with an appropriate chip file.\n\nIf this was intentional, to bypass this check you can set 'override gene list length validation' (--ogllv) to 'True' but this is not recommended."))
    if len(input_ds) < 10000 and options.override == False and options.collapse != "none":
        sys.exit(print("Only ", len(input_ds), "genes were identified in the dataset.\nEither the dataset did not contain all expressed genes, or there was possibly a problem with the chip selected for collapse dataset.\n\nIf this was intentional, to bypass this check you can set 'override gene list length validation' (--ogllv) to 'True' but this is not recommended."))
    if len(input_ds) < 10000 and options.override == True:
        print("Only", len(input_ds), "genes were identified in the dataset, but the user specified overriding this check. Continuing analysis, as-is however this is not recommended. The input dataset should include all expressed genes.")

    # # Parse CLS file
    # labels, phenotypes = GSEAlib.read_cls(options.cls)
    # phenotypes = GSEAlib.match_phenotypes(input_ds, phenotypes)
    # if options.reverse == True and phenotypes.columns[0] == "Labels":
    #     phenotypes["Phenotypes"] = numpy.where((phenotypes["Phenotypes"] == 0) | (
    #         phenotypes["Phenotypes"] == 1), phenotypes["Phenotypes"] ^ 1, phenotypes["Phenotypes"])
    #     labels = {0: labels[1], 1: labels[0]}
    # phenotypes = phenotypes.sort_values('Phenotypes')
    #
    # # Order the dataset using the phenotypes and write out both files
    # input_ds = input_ds.reindex(columns=phenotypes.index)
    input_ds.to_csv('input/gene_by_sample.tsv', sep="\t")
    # pandas.DataFrame(phenotypes['Phenotypes']).transpose().to_csv(
    #     'input/target_by_sample.tsv', sep="\t", index=False)

    # Parse GMT/GMX gene sets files from a list of inputs and create a name:members dict written out as a json file
    if options.gsdb != None:
        with open(options.gsdb) as f:
            gene_sets_dbfile_list = f.read().splitlines()

    gs_data = GSEAlib.read_sets(gene_sets_dbfile_list)
    genesets = gs_data['genesets']
    genesets_descr = gs_data['descriptions']
    with open('input/raw_set_to_genes.json', 'w') as path:
        json.dump(genesets, path,  indent=2)

    # Filter gene sets to just genes in input dataset
    gs_data_subset = GSEAlib.filter_sets(genesets, input_ds.index)
    gs_data_subset_sets = gs_data_subset['genesets']
    gs_data_subset_lengths = gs_data_subset['lengths']
    passing_lengths = dict((key, value) for key, value in gs_data_subset_lengths.items(
    ) if (value >= max(options.min, 1) and value <= options.max))
    passing_sets = {key: gs_data_subset_sets[key]
                    for key in passing_lengths.keys()}
    with open('input/filtered_set_to_genes.json', 'w') as path:
        json.dump(genesets, path,  indent=2)

    # Construct GSEA Settings json file
    gsea_settings = {
        # "number_of_permutations": options.nperm,
        # "permutation": options.perm,
        # "feature_name": "Features",
        # "metric": options.rank_metric,
        # "score_name": options.rank_metric,
        "algorithm": options.method,
        "exponent": options.exponent,
        "maximum_gene_set_size": options.max,
        "minimum_gene_set_size": options.min,
        "remove_gene_set_genes": True,
        # "random_seed": options.seed,
        # "high_text" : str(labels[0]),
        # "low_text" : str(labels[1]),
        "number_of_jobs": options.cpu,
        "number_of_sets_to_plot": options.nplot,
        "gene_sets_to_plot": []
    }

    with open('input/gsea_settings.json', 'w') as path:
        json.dump(gsea_settings, path,  indent=2)

    # Run GSEA
    subprocess.check_output(['gsea', 'data-rank', 'input/gsea_settings.json',
                             'input/gene_by_sample.tsv', 'input/filtered_set_to_genes.json', os.getcwd()])

    # Not Processing Results into figures for ssGSEA (yet?)

    # Zip up results
    if options.zip == True:
        gsea_files = []
        for folderName, subfolders, filenames in os.walk(os.path.relpath(os.getcwd())):
            for filename in filenames:
                # create complete filepath of file in directory
                filePath = os.path.join(folderName, filename)
                # Add file to zip
                gsea_files.append(filePath)
        with ZipFile("gsea_results.zip", "w") as gsea_zip:
            for filename in gsea_files:
                gsea_zip.write(filename)


if __name__ == '__main__':
    main()
