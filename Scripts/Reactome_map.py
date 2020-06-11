#Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import scipy

#Settings
#Species is set because we set the "Project to human" setting in Reactome
species = 'Homo sapiens'
alpha = 0.05 #max p value for quantification

#read the reactome ID to pathway mapping
term_df = pd.read_csv('ReactomePathways.txt', index_col = 0, sep = '\t', header = None)
term_df.columns = ['pathway', 'species']
#subset by species
term_df = term_df[term_df.species == species]

#Load the reactome network
rela = pd.read_csv('ReactomePathwaysRelation.txt',  sep = '\t', header = None)
rela.columns = ['from', 'to']


# Function to find the end term (most broard term)
def find_lowest(id):
    try:
        low = rela[rela['to'] == id]['from'].values[0]
        if low not in rela['from'].values:
            return low
        else:
            return find_lowest(low)
    except:
        return id

# Make mapping between all nodes and their most broard term (takes a while to run)
end_hash = {id:find_lowest(id) for id in list(set(rela['to'].values) | set(rela['from'].values))}


#Map between ID and pathway descriptor
def map_id(id):
    return term_df.loc[id]['pathway']

#Function to count number of significant links per broard category over a set of reactome output files
def run_per_clust(file, clust, alpha):

    react = pd.read_csv(file.format(clust), index_col = 0)
    react = react[react['Entities pValue'].values <= alpha]
    react['end_id'] = [end_hash[id] for id in react.index]

    react['end_path'] = [map_id(id) for id in react.end_id]
    counts = pd.DataFrame(react.groupby('end_path').count()['end_id'])
    counts.columns = ['cluster{}'.format(clust)]
    return counts

#Quantify over 5 clusters
N_nodes = pd.concat([run_per_clust('Cluster{}_result.csv', clust, alpha) for clust in np.arange(1,6,1)], axis = 1).replace(np.nan, 0)

import seaborn as sns
melt = pd.melt(N_nodes.reset_index(), id_vars='index')
melt.columns = ['x', 'y', 'value']


## Function to make a dotplot / Heatmap
# Source: https://www.kaggle.com/drazen/heatmap-with-sized-markers
#           and https://towardsdatascience.com/better-heatmaps-and-correlation-matrix-plots-in-python-41445d0f2bec
def heatmap(x, y, **kwargs):
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = [1]*len(x)

    if 'palette' in kwargs:
        palette = kwargs['palette']
        n_colors = len(palette)
    else:
        n_colors = 256 # Use 256 colors for the diverging color palette
        palette = sns.color_palette("Blues", n_colors)

    if 'color_range' in kwargs:
        color_min, color_max = kwargs['color_range']
    else:
        color_min, color_max = min(color), max(color) # Range of values that will be mapped to the palette, i.e. min and max possible correlation

    def value_to_color(val):
        if color_min == color_max:
            return palette[-1]
        else:
            val_position = float((val - color_min)) / (color_max - color_min) # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            ind = int(val_position * (n_colors - 1)) # target index in the color palette
            return palette[ind]

    if 'size' in kwargs:
        size = kwargs['size']
    else:
        size = [1]*len(x)

    if 'size_range' in kwargs:
        size_min, size_max = kwargs['size_range'][0], kwargs['size_range'][1]
    else:
        size_min, size_max = min(size), max(size)

    size_scale = kwargs.get('size_scale', 500)

    def value_to_size(val):
        if size_min == size_max:
            return 1 * size_scale
        else:
            val_position = (val - size_min) * 0.99 / (size_max - size_min) + 0.01 # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            return val_position * size_scale
    if 'x_order' in kwargs:
        x_names = [t for t in kwargs['x_order']]
    else:
        x_names = [t for t in sorted(set([v for v in x]))]
    x_to_num = {p[1]:p[0] for p in enumerate(x_names)}

    if 'y_order' in kwargs:
        y_names = [t for t in kwargs['y_order']]
    else:
        y_names = [t for t in sorted(set([v for v in y]))]
    y_to_num = {p[1]:p[0] for p in enumerate(y_names)}

    plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1) # Setup a 1x10 grid
    ax = plt.subplot(plot_grid[:,:-1]) # Use the left 14/15ths of the grid for the main plot

    marker = kwargs.get('marker', 's')

    kwargs_pass_on = {k:v for k,v in kwargs.items() if k not in [
         'color', 'palette', 'color_range', 'size', 'size_range', 'size_scale', 'marker', 'x_order', 'y_order'
    ]}

    ax.scatter(
        x=[x_to_num[v] for v in x],
        y=[y_to_num[v] for v in y],
        marker=marker,
        s=[value_to_size(v) for v in size],
        c=[value_to_color(v) for v in color],
        **kwargs_pass_on
    )
    ax.set_xticks([v for k,v in x_to_num.items()])
    ax.set_xticklabels([k for k in x_to_num], #rotation=90,
    horizontalalignment='center',
    fontdict={'fontsize':5})
    ax.set_yticks([v for k,v in y_to_num.items()])
    ax.set_yticklabels([k for k in y_to_num], fontdict={'fontsize':5})

    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    ax.set_facecolor('#F1F1F1')

    # Add color legend on the right side of the plot
    if color_min < color_max:
        ax = plt.subplot(plot_grid[:,-1]) # Use the rightmost column of the plot

        col_x = [0]*len(palette) # Fixed x coordinate for the bars
        bar_y=np.linspace(color_min, color_max, n_colors) # y coordinates for each of the n_colors bars

        bar_height = bar_y[1] - bar_y[0]
        ax.barh(
            y=bar_y,
            width=[5]*len(palette), # Make bars 5 units wide
            left=col_x, # Make bars start at 0
            height=bar_height,
            color=palette,
            linewidth=0
        )
        ax.set_xlim(1, 2) # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
        ax.grid(False) # Hide grid
        ax.set_facecolor('white') # Make background white
        ax.set_xticks([]) # Remove horizontal ticks
        ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3)) # Show vertical ticks for min, middle and max
        ax.yaxis.tick_right() # Show vertical ticks on the right
        plt.subplots_adjust(left=0.40, bottom=0.25, right=None, top=None, wspace=None, hspace=None)

#Make dotplot / heatmap of the quantification
heatmap(
    melt['y'], melt['x'],
    color=melt['value'],
    palette=sns.color_palette("Greys", 256) ,
    size=melt['value'].abs(),
    marker='$\u2713$',
    x_order=N_nodes.columns,
    y_order=sorted(N_nodes.index),
    size_scale=100)
plt.savefig('Reactome_count.pdf')
plt.subtitle('Reactome Therm Counts')
