import os, sys, pandas, numpy

# Simple GMT/GMX to Dict parser
def read_sets(gene_sets_dbfile_list):
    genesets={}
    genesets_descr={}
    for gsdb in gene_sets_dbfile_list:
        gsdb_split = gsdb.split(".")
        if gsdb_split[-1] == "gmt":
            with open(gsdb) as f:
                temp=f.read().splitlines()
                for i in range(len(temp)):
                    gs_line=temp[i].split("\t")
                    gene_set_name=gs_line[0]
                    gene_set_desc=gs_line[1]
                    gene_set_tags=gs_line[2:len(gs_line)]
                    genesets[gene_set_name]=list(set(gene_set_tags))
                    genesets_descr[gene_set_name]=gene_set_desc # Not used yet but should end up in reports eventually
        else:  # is a gmx formatted file
            df_temp=pandas.read_csv(
                gsdb, sep='\t', skip_blank_lines=True).transpose().dropna(how='all')
            for i in range(len(df_temp)):
                gs_line=df_temp.iloc[i][~df_temp.iloc[i].isnull()]
                gene_set_name=gs_line.name
                gene_set_desc=gs_line[0]
                gene_set_tags=gs_line[1:len(gs_line)]
                genesets[gene_set_name]=list(set(gene_set_tags))
                genesets_descr[gene_set_name]=gene_set_desc # Not used yet but should end up in reports eventually
    return genesets


# Simple implementation of a GCT parser
# Accepts a GCT file and returns a Pandas Dataframe with a single index
def read_gct(gct):
    dataset=pandas.read_csv(gct, sep='\t', header=2, index_col=[
        0, 1], skip_blank_lines=True)
    dataset.index.names=["NAME", "Description"]
    dataset_descriptions=dataset.index.to_frame(index=False)
    dataset_descriptions.set_index(["NAME"], inplace=True)
    dataset.index=dataset.index.droplevel(1)  # Drop gene descriptions
    return {'data': dataset, 'row_descriptions': dataset_descriptions["Description"].values}


# Simple implementation of a CHIP Parser for use with ssGSEA
# Reads in a CHIP formatted file and returns a pandas dataframe containing
# the probe to gene mappings
def read_chip(chip):
    chip_df=pandas.read_csv(chip, sep='\t', index_col=0, skip_blank_lines=True)
    return chip_df


# Simple implementation of GSEA DEsktop's Collapse Dataset functions for use
# with ssSGEA
# Accepts an expression dataset in GCT format, a CHIP file, and a
# collapse metric and returns a pandas dataframe formatted version of the
# dataset collapsed from probe level to gene level using the specified metric.
def collapse_dataset(dataset, chip, method="sum"):
    import pandas as pd
    if isinstance(dataset, pandas.DataFrame):
        dataset=dataset
    elif isinstance(dataset, dict) == False:
        dataset=read_gct(dataset)
    if isinstance(chip, pandas.DataFrame) == False:
        chip=read_chip(chip)
    if isinstance(dataset, dict) == True:
        dataset=dataset['data']
    joined_df=chip.join(dataset, how='inner')
    joined_df.reset_index(drop=True, inplace=True)
    annotations=joined_df[["Gene Symbol",
                             "Gene Title"]].drop_duplicates().copy()
    joined_df.drop("Gene Title", axis=1, inplace=True)
    if method.lower() == "sum":
        collapsed_df=joined_df.groupby(["Gene Symbol"]).sum()
    if method.lower() == "mean":
        collapsed_df=joined_df.groupby(["Gene Symbol"]).mean()
    if method.lower() == "median":
        collapsed_df=joined_df.groupby(["Gene Symbol"]).median()
    if method.lower() == "max":
        collapsed_df=joined_df.groupby(["Gene Symbol"]).max()
    if method.lower() == "absmax":
        collapsed_df = joined_df.loc[joined_df.groupby(['Gene Symbol']).idxmax()]
    collapsed_df.index.name="Name"
    return {'data': collapsed_df, 'row_descriptions': annotations["Gene Title"].values}


# Save a GCT result to a file, ensuring the filename has the extension .gct
def write_gct(gct, file_name, check_file_extension=True):
    if check_file_extension:
        file_name = check_extension(file_name, ".gct")

    rows = str(len(gct['data']))
    columns = str(len(gct['data'].columns))

    if len(gct['row_descriptions'])!=int(rows):
        sys.exit("Number of row descriptions ("+ len(gct['row_descriptions'])+ ") not equal to number of row names ("+ rows+ ").")

    row_descriptions = gct['row_descriptions']
    if row_descriptions == None:
        row_descriptions = ['NA']*int(rows)

    m = gct['data'].copy()
    m.insert(loc=0, column='Description', value=gct['row_descriptions'])
    with open(file_name, 'w') as file:
        file.write('#1.2\n' + rows + '\t' + columns + '\n')
        m.to_csv(file, sep='\t', index_label="NAME", mode='w+')
    return(file_name)


# extension e.g. '.gct'
def check_extension(file_name, extension):
    import re
    ext = re.search(extension+"$", file_name.lower())
    if ext == None:
        file_name = file_name + extension
    return file_name


# Read CLS function adapted from https://github.com/broadinstitute/gsea_python/blob/ccal-refactor/gsea/Utils.py
def read_cls(path):
    """
    Returns a Pandas Series with phenotype labels from a .cls
        path (str): Local filepath for a .cls file
        returns (pd.Series): Pandas Series with phenotype labels
    """
    lines = open(path).readlines()
    labels = {
        label: i for i, label in enumerate(lines[1][1:-1].split())
    }
    try:
        labs = lines[2][:-1].split()
        phens = [labels[lab] for lab in lines[2].strip('\n').split()]
        return pandas.concat([pandas.Series(labs, name='Labels'), pandas.Series(phens, name='Phenotypes')], axis=1)
    except KeyError: ## Assume phenotype row is already ints
        phens = list(map(int, lines[2].strip('\n').split()))
        return pandas.concat([pandas.Series(phens, name='Labels'), pandas.Series(phens, name='Phenotypes')], axis=1)


# Read Maps Phenotypes to Samples, adapted from https://github.com/broadinstitute/gsea_python/blob/ccal-refactor/gsea/Utils.py
def match_phenotypes(expr, phen):
    """
    Populates the index of phen with the column names of expr
        expr (pandas.DataFrame): DataFrame from read_gct
        phen (pandas.Series): Series from read_cls
        returns (pandas.Series): phen with index set to expr columns
    """
    if len(phen) == len(expr.columns):
        phen.index = expr.columns
    else:
#        common = set(phen['phenotypes'].index) & set(expr.columns)
#        expr = expr[list(common)]
#        phen = phen[list(common)]
        sys.exit("The number of samples in the CLS file did not match the number of samples in the dataset.")
    return phen
