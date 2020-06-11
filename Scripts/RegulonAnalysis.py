from matplotlib import pyplot as plt
from sklearn import preprocessing
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.colors import to_rgb
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns
import networkx as nx
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.special import rel_entr
from math import sqrt
from math import ceil, floor
try:
    from scipy.spatial.distance import jensenshannon
except:
    def jensenshannon(p, q, base=None):
        """
        Compute the Jensen-Shannon distance (metric) between
        two 1-D probability arrays. 
        """
        p = np.asarray(p)
        q = np.asarray(q)
        p = p / np.sum(p, axis=0)
        q = q / np.sum(q, axis=0)
        m = (p + q) / 2.0
        left = rel_entr(p, m)
        right = rel_entr(q, m)
        js = np.sum(left, axis=0) + np.sum(right, axis=0)
        if base is not None:
            js /= np.log(base)
        return np.sqrt(js / 2.0)



class RegulonsAnalysis():
    """Module for analysis of regulon activity. 
        regulons: Regulon object - Output of RCisTarget.
        scanpy: A scanpy object with gene expression data (must contain 'celltype', 'protocol', 'tissue', 'phase' columns in obs)
        auc_mtx: Output of AUCell
        """
    
    def __init__(self, regulons, scanpy, auc_mtx):
        """Prepare object"""

        self.regulons = regulons
        self.scanpy = scanpy
        self.ex_matrix = pd.DataFrame(
                                data = scanpy.X,
                                index = self.scanpy.obs_names.tolist(),
                                columns = self.scanpy.var.index.tolist())
        
        self.auc_mtx = auc_mtx
        self.remove_regulons()
        self.regulon_scanpy = self.make_regulon_scanpy()
        self.dr = None
        self.pcs = None
    
    def remove_regulons(self):
        """Remove not used regulons"""
        for i in self.auc_mtx.columns:
            if sum(self.auc_mtx[i]) == 0:
                del self.auc_mtx[i]
    

    ######
    def cell_regulon_cor(self, color_index, color_by):
        """Plots spearmman cell-cell correlation"""
        #Calculate spearman correlation of regulon activities of single cells
        self.spear_corr, _ = spearmanr(self.auc_mtx, axis = 1)

        row_colors = self.get_cor_colors(color_index, color_by)
        row_colors = pd.DataFrame(data = row_colors.values, columns=['Celltype'], index = self.scanpy.obs_names.tolist())
        row_colors = pd.DataFrame(data = row_colors.values, columns=['Celltype'], index = self.scanpy.obs_names.tolist())
        if len(set(self.regulon_scanpy.obs.protocol)) > 1: 

            row_colors['protocol'] = self.get_cor_colors_protocol().tolist()


        g = sns.clustermap(pd.DataFrame(self.spear_corr, index = self.scanpy.obs_names.tolist(), columns = self.scanpy.obs_names.tolist()),
                           row_colors=row_colors, 
                           cmap = "viridis", 
                           linewidths=0, 
                           xticklabels=False, 
                           yticklabels=False,
                           vmin= 0, vmax = 1,
                           rasterized = True
                          )
    
    def get_cor_colors(self, color_index, color_by):
        """Helper to get color map"""
        color_map = {cell:color for cell,color in zip(color_index['celltype'], color_index[color_by])}
        color = self.regulon_scanpy.obs['celltype'].map(color_map)
        return color

    def get_cor_colors_protocol(self):

        lut = {'Smartseq2': '#279e68', '10X': '#1f77b4', 'Microwell-seq': '#ff7f0e'}
        color = self.regulon_scanpy.obs['protocol'].map(lut)

        return color

    #######
    def cell_regulon_network(self, thresh, color_index, dr, pcs, color_by, save = None):
        """plot links between highly correlated cells """
        links_filtered = self.link_params(thresh)
        G = self.build_graph(links_filtered)
        posis = self.pos_from_dr(dr, pcs)
        self.plot_network(posis, G, color_index, color_by, save)
        
    
    def link_params(self, thresh):
        """Filter links on correlation threshold"""

        self.spear_corr = pd.DataFrame(data=self.spear_corr, index=self.ex_matrix.index, columns=self.ex_matrix.index)


        links = self.spear_corr.stack().reset_index()
        links.columns = ['var1', 'var2','value']

        # Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
        links_filtered=links.loc[ (links['value'] > thresh) & (links['var1'] != links['var2']) ]
        return links_filtered
 

    def build_graph(self, links_filtered):
        """Build networkx graph"""
        G = nx.from_pandas_edgelist(links_filtered, 'var1', 'var2')
        self.G = G
        return G
    
    def pos_from_dr(self, dr, pcs):
        """Get precomputed embedding """
        #Node position
        posis = dict()
        labels = list()
        dim = self.get_dr(dr, pcs)
        for i in range(len(self.ex_matrix.index)):
            posis[self.ex_matrix.index[i]] = [dim[i,0] , dim[i,1]]
            #labels[self.ex_matrix.index[i]] = self.scanpy.obs['celltype'][i]
            
        return posis
    
    def plot_network(self, posis, G, color_index, color_by, save=None):
        """Plots the network on precomputed embedding and saves network as gexf file"""

        pos=posis
        self.got_colors = False
        # find node near center (0.5,0.5)
        dmin=1
        ncenter=0
        for n in pos:
            x,y=pos[n]
            d=(x-0.5)**2+(y-0.5)**2
            if d<dmin:
                ncenter=n
                dmin=d
                
        for idx, node in enumerate(G.nodes):
            if idx != 0:
                color = self.get_color(G, color_index, color_by)[idx]
                self.got_colors = True
                
            color = self.get_color(G, color_index, color_by)[idx]
            color = to_rgb(color)
            c = []
            for i in color:
                c.append(i)
            c=[int(255*x) for x in c]
            G.nodes[node]['viz'] = {'color': {'r': c[0], 'g': c[1], 'b': c[2], 'a': 1}, 'position':{'x':pos[node][0],'y':pos[node][1]}, 'z':0.0000000000}
            

        plt.figure(figsize=(8,8))
        nx.draw_networkx_edges(G,pos,nodelist=[ncenter],alpha=0.005)
        nx.draw_networkx_nodes(G,pos,
                               node_size=15,
                               node_color=self.get_color(G, color_index, color_by),
                               cmap=plt.cm.Reds_r)
        if type(save) == str:
            nx.write_gexf(G, save)
        
    def get_color(self, G, color_index, color_by):
        """Get colors for the network"""
        if self.got_colors == True:
            colors = self.colors
            return colors
        labels = []
        colors=[]
        for node in G:
            for idx , cell in enumerate(self.auc_mtx.index):
                if cell == node:
                    #find celltype of cell
                    labels.append(self.regulon_scanpy.obs['celltype'].values[idx])
                    

            count = 0
            for x in color_index['celltype'].values:
                if x == labels[-1]:
                    colors.append(color_index[color_by].values[count])
                    
                else:
                    count = count + 1
        self.colors = colors
        return colors
    
    ######
    def pear_corr(self):
        """Calculates pearson correlation for use in 
        calculation"""
    
        self.pear_corr_reg = pd.DataFrame(np.corrcoef(self.auc_mtx.T))

    
    ######
    def calc_csi(self):
        """Calculates Connection specificity index (CSI)"""
        self.pear_corr()
        #Calc CSI
        
        corr_mat = self.pear_corr_reg
        csi_mat = np.zeros(corr_mat.shape)
        
        for i in range(corr_mat.shape[0]):
            a = corr_mat.index[i]
            for j in range(corr_mat.shape[1]):
                b = corr_mat.columns[j]
                c = corr_mat.iloc[i][j] - 0.05
                conn_pais_a = set(corr_mat.index[corr_mat.loc[i , :] >= c])
                conn_pais_b = set(corr_mat.index[corr_mat.loc[: , j] >= c])
                conn_pais_ab = len(conn_pais_a.union(conn_pais_b))
                n = corr_mat.shape[0]
                csi = 1 - (conn_pais_ab / n)
                csi_mat[i,j] = csi
                
        
        self.csi_mat_df = pd.DataFrame(data=csi_mat, index=self.auc_mtx.columns.values, columns= self.auc_mtx.columns.values)
        return self.csi_mat_df
    
    def find_modules(self, csi_mat_df,cluster_n):
        """Performes hierachical clustering """
        clustering = AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
                    connectivity=None, linkage='complete', memory=None, n_clusters=cluster_n).fit(csi_mat_df)
        
        
        correlations = csi_mat_df
        correlations_array = np.asarray(correlations)

        row_linkage = hierarchy.linkage(distance.pdist(correlations_array), method='average')
        col_linkage = hierarchy.linkage(distance.pdist(correlations_array.T), method='average')
        clustering = hierarchy.fcluster(row_linkage,cluster_n, 'maxclust')
 
        return clustering
    
    def plot_csi(self, cluster_n = 5, cluster = None, lut = None):
        """Plots the CSI matrix. 
            Performes hierachical clustering if no clustering in input. 
            If inputting clustering, one should also input colormap as lut"""

        csi_mat_df = self.calc_csi()
        
        if cluster is None:
            clustering = self.find_modules(csi_mat_df, cluster_n)
            self.clustering = clustering
            lut = dict(zip(set(clustering), sns.color_palette("husl", len(set(clustering)))))

            self.lut = lut
        else:
            self.clustering = np.array(cluster, dtype='int32')
            self.lut = lut
            
        row_colors = pd.DataFrame(self.clustering)[0].map(self.lut)
        self.cluster_row_colors = row_colors
        
        g = sns.clustermap(csi_mat_df,
                           row_colors=row_colors.values,
                           figsize=(20,20),
                           xticklabels = [],#self.get_new_names(csi_mat_df.index.tolist(), self.regulons),#self.get_new_names(self.get_regulon_names(self.regulons), self.regulons), 
                           yticklabels = [],#self.get_new_names(csi_mat_df.index.tolist(), self.regulons), #self.get_new_names(self.get_regulon_names(self.regulons), self.regulons),
                           cmap ="viridis",  
                           annot_kws={"size": 3},
                          rasterized = True
                          )

        
    #######    
    def regulon_network(self, thresh, save, lut=None):
        """Plots the regulon-regulon network based ont he CSI matrix"""
        self.calc_links(thresh)
        G = self.build_graph_reg()
        self.plot_reg_network(G, self.color_nodes(G, lut), save)
        
        
    def calc_links(self, thresh):
        """Subsets link matrix based on threshold"""
        
        if hasattr(self, 'csi_mat_df'):
            pass
        else:
            self.calc_csi()
        
        #regulon network
        self.links_reg = self.csi_mat_df.stack().reset_index()
        self.links_reg.columns = ['var1', 'var2','value']

        # Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
        self.links_filtered_reg=self.links_reg.loc[ (self.links_reg['value'] > thresh) & (self.links_reg['var1'] != self.links_reg['var2']) ]
    
    
    def build_graph_reg(self):
        """Build networks graph"""
        G=nx.from_pandas_edgelist(self.links_filtered_reg, 'var1', 'var2')
        return G
        
    def color_nodes(self, G, lut):
        """returns color labels of regulons"""
        labels = []
        
        if not hasattr(self, 'clustering'):
            self.clustering = self.find_modules(self.calc_csi())
        
        for idx, node in enumerate(G):
            labels.append(self.clustering[self.auc_mtx.columns.values == node])
            G.nodes[node]['cluster'] = str(self.clustering[self.auc_mtx.columns.values == node])
            
            self.col = self.cluster_row_colors[self.auc_mtx.columns.values == node]
            c = []
            for i in self.col:
                c.append(i[0])
                c.append(i[1])
                c.append(i[2])
            c=[int(255*x) for x in c]
            G.nodes[node]['viz'] = {'color': {'r': c[0], 'g': c[1], 'b': c[2], 'a': 1}}

        if lut is None:    
            lut = dict(zip(set(self.clustering), sns.color_palette("husl", len(set(self.clustering)))))
        labels = pd.DataFrame(labels)[0].map(lut)
        
        return labels
        
        
    def plot_reg_network(self, G, labels, save ):
        """Plots regulon network and stores network as gexf file"""
        nx.draw(G, 
                with_labels=True, 
                node_size=200,
                node_color=labels, 
                alpha=0.5)
        if type(save) == str:
            nx.write_gexf(G, save)
            
    ######
    def avg_reg_activity(self, celltypes):
        #calculate average regulon activity of a celltype
        avg_activity = pd.DataFrame()
        
        if celltypes != ['all']:
            types = celltypes
        else:
            types = set(self.scanpy.obs['celltype'])
            
        for celltypes in types:
            a = []
            #iterates over unique celltypes
            for cell in self.scanpy.obs['celltype']:
                #finding all cells of the given celltype
                a.append(cell == celltypes)
                #Calculate enrichmean as the mean over all cells of the celltype
            avg_activity[celltypes] = self.auc_mtx[a].mean()
            #print(avg_activity.T[])

        for i in avg_activity.index:
            if (sum(avg_activity.T[i]) == 0) or (avg_activity.T[i].std == 0):

                avg_activity = avg_activity.drop(i, axis=0)

        return avg_activity
        
    def plot_avg_reg_activity(self, z_score=None, scale=None,  row_cluster = True):

            
        avg_reg = self.avg_reg_activity(celltypes)
        g = sns.clustermap(avg_reg.T, z_score=z_score,row_cluster=row_cluster, col_cluster=True, standard_scale=scale, cmap ="viridis", figsize=(20,20),
                          rasterized = True)#, xticklabels = self.get_new_names(self.get_regulon_names(self.regulons), self.regulons))
        #g.savefig(self.output_dir + save , dpi=1200, transparent=True)
    #######    
    def binarize(self):
        from pyscenic.binarization import binarize
        from pyscenic.binarization import plot_binarization
        #Binarize AUCell_mtx
        binar, auc_thresholds = binarize(self.auc_mtx)

        #plot_binarization(auc_mtx=self.auc_mtx, regulon_name="Cebpb(+)", bins = 25, threshold = None)#auc_thresholds)
        return binar
    
    def reg_enrichment(self):
        #Calculate enrichment of regulon activity in each celltype
        binar = self.binarize()
        reg_enrichment = pd.DataFrame()
        for celltypes in set(self.scanpy.obs['celltype']):
            a = []
            #iterates over unique celltypes
            for cell in self.scanpy.obs['celltype']:
                #finding all cells of the given celltype
                a.append(cell == celltypes)
            #Calculate enrichmean as the mean over all cells of the celltype
            reg_enrichment[celltypes] = (binar[a].sum()) / len(binar[a])
        return reg_enrichment
    
    def plot_reg_enrichment(self, z_score=None, scale=None):
        reg_enrichment = self.reg_enrichment()
        sns.clustermap(reg_enrichment.T, z_score=z_score,row_cluster=True, cmap ="viridis", figsize=(20,20), xticklabels = self.get_new_names(self.get_regulon_names(self.regulons), self.regulons) , rasterized = True)
        
        
        
    #######    
    def make_regulon_scanpy(self):
        """Make scanpy object with on AUCell regulon activity matrix"""
        self.auc_mtx.index = self.scanpy.obs_names.tolist()
        test = sc.AnnData(X = self.auc_mtx , obs = self.scanpy.obs['celltype'])
        test.obs['celltype'] =self.scanpy.obs['celltype'].values
        test.obs['protocol'] =self.scanpy.obs['protocol'].values
        test.obs['tissue'] =self.scanpy.obs['tissue'].values
        test.obs['phase'] =self.scanpy.obs['phase'].values
        try:
            test.obs['Super'] =self.scanpy.obs['Super'].values
        except:
            pass
        
        if 'original_annotation' in self.scanpy.obs.columns:
            test.obs['original_annotation'] =self.scanpy.obs['original_annotation'].values
        return test
        
    #######    
    def new_pca(self, plot, color = 'celltype'):
        if color == None:
            color = 'celltype'
        self.regulon_scanpy = self.make_regulon_scanpy()
        sc.tl.pca(self.regulon_scanpy, n_comps=20)
        if plot == True:
            sc.pl.pca(self.regulon_scanpy, color = [color] )
        else:
            pass

    
    ######
    def color_by_celltype(self, color_index, color_by):
        """Make color map to use cell  color scheme """
        colors = {cell:color for cell,color in zip(color_index['celltype'], color_index[color_by])}

        return colors
    
    

    
    def get_dr(self, dr,pcs):
        """Computes and returns embedding"""
        if dr != self.dr and pcs != self.pcs:
            self.dr = dr
            self.pcs = pcs  
            if dr == 'pca':
                self.comp_pca(pcs)
                return self.regulon_scanpy.obsm['X_pca']
            if dr == 'tsne':
                self.comp_tsne(pcs)
                return self.regulon_scanpy.obsm['X_tsne']
            if dr == 'umap':
                if 'X_pca' in  self.regulon_scanpy.obsm:
                    self.comp_umap(pcs)
                    return self.regulon_scanpy.obsm['X_umap']
                else:
                    print('PCA not computed - using %s principal components' % pcs)
                    self.comp_pca(pcs)
                    self.comp_umap(pcs)
                    return self.regulon_scanpy.obsm['X_umap']
        else:
            return self.regulon_scanpy.obsm['X_'+dr]
            
    def comp_pca(self,pcs):
        sc.tl.pca(self.regulon_scanpy, n_comps=pcs)
        
    def comp_tsne(self,pcs):
        sc.tl.tsne(self.regulon_scanpy,n_pcs=pcs)
        
    def comp_umap(self,pcs):
        sc.pp.neighbors(self.regulon_scanpy, n_neighbors=10)
        sc.tl.umap(self.regulon_scanpy)

        
    
    def plot_by_celltype(self, dr, color_index, pcs, color_by, point_size, Type = None, title = None):
        """Plot embedding and color by celltype"""

        dr_mat = self.get_dr(dr,pcs)
        colors = self.color_by_celltype(color_index, color_by)
        plotting_df = pd.DataFrame()
        
        plotting_df['X_dr'] = list(dr_mat[:,0])
        plotting_df['Y_dr'] = list(dr_mat[:,1])
        plotting_df['type'] = pd.Categorical(self.regulon_scanpy.obs['celltype'])
        sns.set_style( "white", {'axes.grid': False})
        plt.scatter(
            x=plotting_df['X_dr'],
            y = plotting_df['Y_dr'], 
            c = self.regulon_scanpy.obs['celltype'].map(colors),
            s=point_size)
        plt.gca().axes.yaxis.set_ticklabels([])
        plt.gca().axes.xaxis.set_ticklabels([])
        plt.xlabel(dr.upper() + '1')
        plt.ylabel(dr.upper() + '2')
        if title == None:
            title = 'celltype'
        plt.title(title)

    
    def plot_by_celltype_ind(self, Type, dr, color_index,pcs, point_size=100):
        #transform_color_index
        bo = color_index['celltype'] == Type
        color_index['ind'] = [d if t == True else '#E0E0E0' for t, d in zip(bo, color_index['color_grad'])]
        
        self.plot_by_celltype(dr, color_index,pcs, color_by='ind', point_size = point_size, Type = Type, title = Type)
        #should drop ind column after running
    
    def plot_activity(self, regulon, dr, pcs):
        """Plots regulon activity on embedding"""
        self.get_dr(dr, pcs)
        if dr == 'pca':
            sc.pl.pca(self.regulon_scanpy, color = regulon, color_map="viridis")
        if dr == 'tsne':
            sc.pl.tsne(self.regulon_scanpy, color = regulon, color_map="viridis")
        if dr == 'umap':
            sc.pl.umap(self.regulon_scanpy, color = regulon, color_map="viridis")
        
        
    def plot_by_cellcycle(self,dr, color_index,pcs, color_index_type ):
        """Plot cell cycle phase on embedding"""

        if dr == 'pca':
            sc.pl.pca(self.regulon_scanpy, color = 'phase', palette=['#c4c2c2','#ff4c4c','#f29f9f'])
        if dr == 'tsne':
            sc.pl.tsne(self.regulon_scanpy, color = 'phase', palette=['#c4c2c2','#ff4c4c','#f29f9f'])
        if dr == 'umap':
            sc.pl.umap(self.regulon_scanpy, color = 'phase', palette=['#c4c2c2','#ff4c4c','#f29f9f'])
        
    
    def plot(self, dr, color_by, color_index=None, pcs=5, color_index_type='color', point_size=100):
        if color_by == 'celltype':
            self.plot_by_celltype(dr, color_index,pcs, color_index_type, point_size)
        elif color_by == 'cellcycle':
            self.plot_by_cellcycle(dr, color_index,pcs, color_index_type)
        else:
            self.plot_activity(color_by, dr, pcs)
    
    



    def regulon_specificity_scores(self, feature):
        """
        Calculates the Regulon Specificty Scores (RSS). [doi: 10.1016/j.celrep.2018.10.045]
        :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
        :param cell_type_series: A pandas Series object with cell identifiers as index and cell type labels as values.
        :return: A pandas dataframe with the RSS values (cell type x regulon).

        from https://github.com/aertslab/pySCENIC/blob/master/src/pyscenic/rss.py

        """
        cell_type_series = self.regulon_scanpy.obs[feature]
        cell_types = list(cell_type_series.unique())
        n_types = len(cell_types)
        regulons = list(self.auc_mtx.columns)
        n_regulons = len(regulons)
        rss_values = np.empty(shape=(n_types, n_regulons), dtype=np.float)

        def rss(aucs, labels):
            # jensenshannon function provides distance which is the sqrt of the JS divergence.
            return 1.0 - jensenshannon(aucs/aucs.sum(), labels/labels.sum())

        for cidx, regulon_name in enumerate(regulons):
            for ridx, cell_type in enumerate(cell_types):
                rss_values[ridx, cidx] = rss(self.auc_mtx[regulon_name], (cell_type_series == cell_type).astype(int))

        self.rss = pd.DataFrame(data=rss_values, index=cell_types, columns=regulons)

    def plot_regulon_specificity_celltype(self, cell_type, top_n=5, max_n=None):
        if not hasattr(self, 'rss'):
            self.regulon_specificity_scores('celltype')
        if max_n is None:
            max_n = self.rss.shape[1]

        data = self.rss.T[cell_type].sort_values(ascending=False)[0:max_n]
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
        ax.plot(np.arange(len(data)), data, '.')
        ax.set_ylim([floor(data.min()*100.0)/100.0, ceil(data.max()*100.0)/100.0])
        ax.set_ylabel('RSS')
        ax.set_xlabel('Regulon')
        ax.set_title(cell_type)
        ax.set_xticklabels([])

        font = {
            'color':  'red',
            'weight': 'normal',
            'size':16 ,
        }

        for idx, (regulon_name, rss_val) in enumerate(zip(data[0:top_n].index, data[0:top_n].values)):
            ax.plot([idx, idx], [rss_val, rss_val], 'r.')
            ax.text(idx+(max_n/25), rss_val, regulon_name, fontdict=font, horizontalalignment='left', verticalalignment='center')

        
    def get_new_names(self, regulon_name_list, regulons ):
        a = []
        for i in regulon_name_list:
            a.append(i + ' ' + str(len(self.find_reg_genes(i, regulons))) + 'g')
        return a
    
    def get_regulon_names(self, regulon):
        regulon_names = []
        for i in regulon:
            regulon_names.append(i.name)
        return regulon_names
    
    def find_reg_genes(self, regulon_name, regulons):
        x=0
        for i in regulons:

            if i.name != regulon_name:
                x = x + 1
            else:
                return [k  for  k in  regulons[x].gene2weight.keys()]