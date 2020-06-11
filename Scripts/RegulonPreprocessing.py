import pandas as pd
import scanpy.api as sc
import scanpy
import numpy as np
def tissue_from_name(scanpy):
    """Extracts tissue information from pseudobulk cell names"""
    
    tissue = list()
    for idx, i  in enumerate(scanpy.obs.index):
        if i.startswith('Fetal_'):
            tissue.append((i.split('_')[0])+'_'+(i.split('_')[1]))
        elif i.startswith('Large_'):
            tissue.append((i.split('_')[0])+'_'+(i.split('_')[1]))
        else:
            tissue.append(i.split('_')[0])
        #also put in tissue
        #print(len(tissue))
    scanpy.obs['tissue'] = tissue
    scanpy.obs['tissue'] = pd.Categorical(scanpy.obs['tissue'])

#score cell cycle
def convert_genes(l):
    """Converts regev_lab_cell_cycle_genes.txt to first char uppercase, and the rest lowercase"""
    new = []
    for i in l:
        new.append(i[0]+i[1:].lower())
    return new

def score_cell_cycle(scanpy):

    """Runs the scanpy score_genes_cell_cycle method using the regev_lab_cell_cycle_genes.txt gene list"""
    cell_cycle_genes = [x.strip() for x in open('data/regev_lab_cell_cycle_genes.txt')]
    cell_cycle_genes = convert_genes(cell_cycle_genes)

    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in scanpy.var_names]

    sc.tl.score_genes_cell_cycle(scanpy, s_genes=s_genes, g2m_genes=g2m_genes)

    adata_cc_genes = scanpy[:, cell_cycle_genes]
    sc.tl.pca(adata_cc_genes)

def load_PB_datasets(path = 'data'):
    #Import datasets
    '/'.join([path, 'droplet_pseudobulk.h5ad'])
    droplet = scanpy.read_h5ad('/'.join([path, 'droplet_pseudobulk.h5ad']))
    facs = scanpy.read_h5ad('/'.join([path, 'facs_pseudobulk.h5ad']))
    mca = scanpy.read_h5ad('/'.join([path, 'mca_pseudobulk.h5ad']))
    return droplet, facs, mca

def retain_overlapping_genes(droplet, facs, mca):
    #Only use overlap between all datasets
    overlap = set(droplet.var_names.values) & set(facs.var_names.values) & set(mca.var_names.values)


    droplet = droplet[:, list(overlap)]
    facs = facs[:, list(overlap)]
    mca = mca[:, list(overlap)]

    sc.pp.filter_genes(droplet, min_cells=(len(droplet.obs_names) / 10))
    sc.pp.filter_genes(facs, min_cells=(len(facs.obs_names) / 10))
    sc.pp.filter_genes(mca, min_cells=(len(mca.obs_names) / 10))

    overlap = set(droplet.var_names.values) & set(facs.var_names.values) & set(mca.var_names.values)


    droplet = droplet[:, list(overlap)]
    facs = facs[:, list(overlap)]
    mca = mca[:, list(overlap)]
    return droplet, facs, mca

def tissue_from_name(scanpy):
    tissue = list()
    for idx, i  in enumerate(scanpy.obs.index):
        if i.startswith('Fetal_'):
            tissue.append((i.split('_')[0])+'_'+(i.split('_')[1]))
        elif i.startswith('Large_'):
            tissue.append((i.split('_')[0])+'_'+(i.split('_')[1]))
        else:
            tissue.append(i.split('_')[0])
    scanpy.obs['tissue'] = tissue
    scanpy.obs['tissue'] = pd.Categorical(scanpy.obs['tissue'])

def preprocess_pb(data_path = 'data'):
    print('Loading datasets')
    droplet, facs, mca = load_PB_datasets(data_path)
    print('Subsetting genes (expressed in 10% of cells in all datasets)')
    droplet, facs, mca = retain_overlapping_genes(droplet, facs, mca)

    #library size normalization
    print('Normalizing')
    sc.pp.normalize_per_cell(droplet, counts_per_cell_after=1e4)
    sc.pp.normalize_per_cell(facs, counts_per_cell_after=1e4)
    sc.pp.normalize_per_cell(mca, counts_per_cell_after=1e4)

    #log transform
    print('Log transforming')
    sc.pp.log1p(droplet)
    sc.pp.log1p(facs)
    sc.pp.log1p(mca)

    #identify highly variable genes droplet
    print('Detecting highly variable genes')
    sc.pp.highly_variable_genes(droplet, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.highly_variable_genes(facs, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.highly_variable_genes(mca, min_mean=0.0125, max_mean=3, min_disp=0.5)

    #find mito genes
    mito_genes_drop = droplet.var_names.str.startswith('Mt')
    mito_genes_facs = facs.var_names.str.startswith('Mt')
    mito_genes_mca = mca.var_names.str.startswith('Mt')

    # for each cell compute fraction of counts in mito genes vs. all genes
    droplet.obs['percent_mito'] = np.sum(
        droplet[:, mito_genes_drop].X) / np.sum(droplet.X)
    facs.obs['percent_mito'] = np.sum(
        facs[:, mito_genes_facs].X) / np.sum(facs.X)
    mca.obs['percent_mito'] = np.sum(
        mca[:, mito_genes_mca].X) / np.sum(mca.X)

    #regress out unwanted sources of variation
    print('Regressing out: # counts and % mito')
    sc.pp.regress_out(droplet, ['n_counts', 'percent_mito'])
    sc.pp.regress_out(facs, ['n_counts', 'percent_mito'])
    sc.pp.regress_out(mca, ['n_counts', 'percent_mito'])   

    #scale to unit variance
    print('Scaling')
    sc.pp.scale(droplet, max_value=10)
    sc.pp.scale(facs, max_value=10)
    sc.pp.scale(mca, max_value=10)

    #remove unassigned cells
    print('Removing unannotated cells')
    droplet = droplet[droplet.obs['celltype'] != 'unassigned']
    facs = facs[facs.obs['celltype'] != 'unassigned']
    mca= mca[mca.obs['celltype'] != 'unassigned']

    #annotate protocol
    print('Annotating datasets')
    droplet.obs['protocol'] = '10X'
    facs.obs['protocol'] = 'Smartseq2'
    mca.obs['protocol'] = 'Microwell-seq'

    tissue_from_name(droplet)
    tissue_from_name(facs)
    tissue_from_name(mca)

    try:
        score_cell_cycle(droplet)
    except:
        pass

    try:
        score_cell_cycle(facs)
    except:
        pass

    try:
        score_cell_cycle(mca)
    except:
        pass

    return droplet, facs, mca


def info_table(regulons):
    regulon_names = [x.name for x in regulons]
    #function that outputs a list of gene names present in a given regulon
    def find_reg_genes(regulon_name, regulons):
        x=0
        for i in regulons:
            if i.name != regulon_name:
                x = x + 1
            else:
                return [k  for  k in  regulons[x].gene2weight.keys()]

    genes = [find_reg_genes(x, regulons) for x in regulon_names]
    regulon_info = pd.DataFrame({'Regulon':[x.strip('(+)') for x in regulon_names],
                  'genes':genes,
                  'N_genes':[len(x) for x in genes] })
    regulon_info.index = regulon_info['Regulon'].astype(str)
    return regulon_info

def preprocess_single_cell():
    print('Loading datasets')
    droplet = sc.read('data/tm_droplet_scanpy_no_processing.h5ad')
    facs = sc.read('data/tm_facs_scanpy_no_processing.h5ad')
    mca = sc.read('data/mca_scanpy_no_processing.h5ad')

    print('Getting annotation')
    mca_ann = pd.read_csv('data/mca_annotation_projected_from_drop.tsv', sep= '\t')
    facs_ann = pd.read_csv('data/facs_annotation_projected_from_drop.tsv', sep= '\t')



    mca.obs['original_annotation'] =  mca_ann['original'].values
    mca.obs['from_droplet'] =  mca_ann['projected_drom_droplet'].values

    facs.obs['original_annotation'] =  facs_ann['original'].values
    facs.obs['from_droplet'] =  facs_ann['projected_drom_droplet'].values

    droplet.obs['protocol'] = '10X'
    facs.obs['protocol'] = 'Smartseq2'
    mca.obs['protocol'] = 'Microwell-seq'

    print('Removing cells not annotated by author')

    droplet = droplet[[type(i) == str for i in droplet.obs['cell_ontology_class']]]
    droplet.obs['cell_ontology_class'].replace(np.NaN, 'Non-annotated', inplace = True)
    facs.obs['cell_ontology_class'].replace(np.NaN, 'Non-annotated', inplace = True)
    mca.obs['original_annotation'].replace(np.NaN, 'Non-annotated', inplace = True)

    droplet = droplet[droplet.obs['cell_ontology_class'] != 'Non-annotated']
    facs = facs[facs.obs['cell_ontology_class'] != 'Non-annotated']
    mca = mca[mca.obs['original_annotation'] != 'Non-annotated']

    droplet.obs['celltype'] = droplet.obs['cell_ontology_class']
    droplet = droplet[droplet.obs['celltype'] != 'unassigned']
    droplet.obs['original_annotation'] = droplet.obs['cell_ontology_class']


    for i in [facs, droplet, mca]:
        i.obs['n_counts'] = i.X.sum(axis=1).A1
        sc.pp.calculate_qc_metrics(i, inplace = True)


    print('Removing cells unassigned by scmap')
    droplet.obs['celltype'] = droplet.obs['cell_ontology_class']
    droplet = droplet[droplet.obs['celltype'] != 'unassigned']

    #droplet = droplet[[type(i) == str for i in droplet.obs['celltype']]]

    facs.obs['celltype'] = pd.Categorical(facs.obs['from_droplet'])
    facs = facs[facs.obs['celltype'] != 'unassigned']
    facs = facs[[type(i) == str for i in facs.obs['cell_ontology_class']]]

    mca.obs['celltype'] = pd.Categorical(mca.obs['from_droplet'])
    mca = mca[mca.obs['celltype'] != 'unassigned']
    mca = mca[[type(i) == str for i in mca.obs['celltype']]]

    print('Classifying cell-cycle phase')

    try:
        score_cell_cycle(droplet)
    except:
        pass

    try:
        score_cell_cycle(facs)
    except:
        pass

    try:
        score_cell_cycle(mca)
    except:
        pass
    
    droplet.obs['n_counts'] = droplet.X.sum(axis=1).A1
    facs.obs['n_counts'] = facs.X.sum(axis=1).A1
    mca.obs['n_counts'] = mca.X.sum(axis=1).A1
    
    

    print('Filtering genes')
    sc.pp.filter_genes(droplet, min_counts=1)
    sc.pp.filter_genes(facs, min_counts=1)
    sc.pp.filter_genes(mca, min_counts=1)
    
    
    #library size normalization
    print('Normalizing to library size')
    sc.pp.normalize_per_cell(droplet, counts_per_cell_after=1e4)
    sc.pp.normalize_per_cell(facs, counts_per_cell_after=1e4)
    sc.pp.normalize_per_cell(mca, counts_per_cell_after=1e4)
    
    print('Log transforming')
    #log transform
    sc.pp.log1p(droplet)
    sc.pp.log1p(facs)
    sc.pp.log1p(mca)

    mito_genes_drop = droplet.var_names.str.startswith('Mt')
    mito_genes_facs = facs.var_names.str.startswith('Mt')
    mito_genes_mca = mca.var_names.str.startswith('Mt')

    droplet.obs['percent_mito'] = np.sum(
        droplet[:, mito_genes_drop].X) / np.sum(droplet.X)
    facs.obs['percent_mito'] = np.sum(
        facs[:, mito_genes_facs].X) / np.sum(facs.X)
    mca.obs['percent_mito'] = np.sum(
        mca[:, mito_genes_mca].X) / np.sum(mca.X)

    print('Regressing out library size and percent mitochondrial counts')
    #Skipped to spead up execution time. This was applied when running Aucell in the loaded precompued auc_mtx.csv files
    #regress out unwanted sources of variation
    #sc.pp.regress_out(droplet, ['n_counts', 'percent_mito'])
    #sc.pp.regress_out(facs, ['n_counts', 'percent_mito'])
    #sc.pp.regress_out(mca, ['n_counts', 'percent_mito'])
    
    print('Scaling genes to mean 0 and unit variance')
    #scale to unit variance
    sc.pp.scale(droplet, max_value=10)
    sc.pp.scale(facs, max_value=10)
    sc.pp.scale(mca, max_value=10)


    return droplet, facs, mca

def stats(obj, name):
    df = pd.DataFrame()
    for idx, i in enumerate(obj):
        df.loc['# Tissues', name[idx]] = len(set(i.obs.tissue)) #n_tissue
        df.loc['Mean Counts', name[idx]] =np.mean(i.obs.n_counts)
        df.loc['Median Counts', name[idx]] =np.median(i.obs.n_counts)
        
        df.loc['Mean Genes', name[idx]] =np.mean(i.obs.n_genes_by_counts)
        df.loc['Median Genes', name[idx]] =np.median(i.obs.n_genes_by_counts)
        
        if 'original_annotation' in i.obs.columns:
            df.loc['Author cell-types', name[idx]] =len(set(i.obs['original_annotation']))
        else:
            df.loc['Author cell-types', name[idx]] =len(set(i.obs['cell_ontology_class']))
        df.loc['Cells', name[idx]] =i.X.shape[0]
    return df
