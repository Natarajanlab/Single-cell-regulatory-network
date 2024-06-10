# Single-cell-regulatory-network

This repository contain code to infer and analyse gene regulatory networks (GRN's) in single cell RNA-seq atlases Tabula Muris and Mouse Cell Atlas using SCENIC. 

## Structure
* Data/Regulons
* Notebooks
* Scripts

The data directory contain the inferred regulons. Raw data is not included due to file size limitaitons. 

The notebooks directory contains all notebooks describing the analysis.
The notebooks are split into 5 groups:
1. Prep_objects (Modifying raw data for import into R and Python enviorments)
2. scmap (Notebooks for mapping celltypes between atlases)
3. Make_pseudobulk (notebook to generate pseudobulk datasets)
4. Regulon_inference (Notebooks describing regulon inference with SCENIC)
5. Regulon_analysis (Notebooks describing the analysis of inferred GRN's)

The scripts directory contains misc scripts and functions used throughout the notebooks.

## Citation
We are pleased to share that this work is now published.
**Predicting gene regulatory networks from cell atlases**
Møller AF, Natarajan KN. 
Life Sci Alliance. 2020 Sep 21;3(11):e202000658. doi: 10.26508/lsa.202000658. PMID: 32958603; PMCID: PMC7536823.
_https://www.life-science-alliance.org/content/3/11/e202000658_
