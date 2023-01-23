## core
from datetime import datetime
import time
import os
import ipdb
import sys
# set directory to where the data lies
# os.chdir("/home/rstudio/data/Ahmad_workdir")

## auxiliary
import scanpy as sc
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import importlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D

## via
import pyVIA.core as via
from pyVIA.examples import via_wrapper
import draw_traj # py script in the same directory
importlib.reload(draw_traj)
from draw_traj import draw_trajectory_gams_f

## palantir
import palantir
import plot_palantir # py script in the same directory
importlib.reload(plot_palantir)
from plot_palantir import plot_palantir_results_f

## paga
import scanpy.api as sc # conflicts with above sc import?
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')  # low dpi (dots per inch) yields small inline figures
# sc.logging.print_versions()
results_file = './write/out.h5ad'

def evaluate_speed(data, labels, dimred, method=['via'],
                   origin = None, ori_id = 0, sizes=[1000, 'all'], repeat = 1,
                   quality='coarse', via_allcurves = False, ncomps = -1):
                     
    """Function used to test one or multiple trajectory inference methods, 
    returning running times, post-analysis objects, and plotting the result.

    Parameters
    ----------
    data : pandas dataframe
        Flow, CyTOF or scRNASeq data in numeric dataframe of shape n cells x m markers.
        The easiest way to create a pandas dataframe is saving the data to .csv then using pandas.read_csv().
    labels : pandas series
        Annotations for each cell of data. 
        Can also be created by saving a label list to .csv then using pandas.read_csv().
    dimred : pandas dataframe
        Low-dimensional embeddings of data.
    method : list of str, default ['via']
        Names of the trajectory inference methods to evaluate. 
        Currently available: 'via', 'palantir', 'paga', 'stream'
    origin : str, optional (default None)
        Label of the cell population to be used as origin to compute trajectories. 
    ori_id : int, optional (default 0)
        Numeric index of the cell to be used as origin to compute trajectories. 
        Specify either origin or ori_id but not both.
    sizes : list, optional (default [1000, 'all'])
        Size of the dataset(s) to feed to the method. 
        Can accept multiple values to test on multiple datasets. 
        Datasets for a value smaller then the original size of data are prepared via uniform random sampling. 
        Origin cell is then added to the dataset.
        If 'all', takes the entire dataset without any sampling. 
        Do not set any value to a multiple of 50,000 (causes a very specific issue in one method)
    repeat : int, optional (default 1)
        Amount of times to repeat the entire process. 
        Set to 5-10 for a more accurate measure of speed (will take longer).
    quality : str, optional (default coarse)
        Specific parameter for via. 
        Controls the quality of the clustering, pseudotime values and lineage probabilities for via.
        Set to 'coarse' or 'fine'.
    via_allcurves : bool, optional (default False)
        Specific parameter for via.
        If True, will plot all computed trajectories on the resulting graph (may clutter).
    ncomps : int, optional (default -1)
        Number of pca dimensions to use for neighbor graphs, diffusion maps, etc.
        Usually set to 30-100.

    Returns
    -------
    tuple
        a tuple with two values:
        times: dict
            Each key corresponds to a method. 
            The associated value is a dict where each key corresponds to a dataset size.
            The values are the average running times for each dataset tested.
        objects: dict
            Each key corresponds to a method. 
            The associated value is a dict where each key corresponds to a dataset size.
            The values are the objects returned by the test functions associated with the tested methods. 
    """
    
    objects = {}
    n = data.shape[0]
    times = dict.fromkeys(sizes)
    
    labels.index = data.index
    dimred.index = data.index

    # sample datasets and find and append origin
    for s in sizes:
        if ori_id is None:
            ori_id = np.where(labels==origin)[0][0]
        if s=='all':
            ss = np.arange(n)
        else:
            ss = np.random.choice(np.arange(n), s, replace=False)
        if ori_id not in ss: 
            ss = np.append(ss, ori_id)
            ori_id_new = ss.shape[0] - 1
        else:
            ori_id_new = np.where(ss==ori_id)[0][0]

        sdata = data.iloc[ss]
        sdimred = dimred.iloc[ss]
        slabels = (labels.iloc[ss]).values.tolist()
        if isinstance(slabels[0], list):
            slabels = [l[0] for l in slabels]
        
        origin_new = ori_id_new

        ###
        # test methods
        ###
        
        st = 0
        for i in range(repeat):
            print("testing {} events, iter {}".format(s, i))
            sys.stdout.flush()
            
            if 'via' in method:
                if not isinstance(sdimred, np.ndarray):
                    sdimred = sdimred.to_numpy()
                print('starting via')
                start = time.perf_counter()
                res = test_via(sdata, slabels, origin_new, quality, ncomps)
                if method in objects:
                    objects[method][s] = res, sdimred
                else:
                    objects[method] = dict.fromkeys(sizes)
                    objects[method][s] = res, sdimred
                    
                end = time.perf_counter()
                draw_trajectory_gams_f(via_coarse=res[0],via_fine=res[0], 
                     embedding=sdimred, draw_all_curves=via_allcurves,
                     scatter_size=10, scatter_alpha=1,)
                
                now = datetime.now().strftime('%m-%d-%H%M%S')
                plt.savefig('./plots/'+method+'_'+str(s)+'_traj_'+'coarse'+'_'+now)
                plt.close('all')    
                
                if quality=='fine':
                    # Finer adjustments to the plot can be made using the arguments to this function
                    draw_trajectory_gams_f(via_coarse=res[0],via_fine=res[1], 
                        embedding=sdimred, draw_all_curves=via_allcurves,
                        scatter_size=10, scatter_alpha=1,)
                    plt.savefig('./plots/'+method+'_'+str(s)+'_traj_'+'fine'+'_'+now)
                    plt.close('all')    

                
                via.via_streamplot(via_coarse=res[0], embedding=sdimred)
                plt.savefig('./plots/'+method+'_'+str(s)+'_streamplot_'+'_'+now)
                plt.close('all')    
                
            elif 'palantir' in method:
                print('starting palantir')
                start = time.perf_counter()
                res = test_palantir(sdata, sdimred, origin_new, n_pca=ncomps)
                if method in objects:
                    objects[method][s] = res
                else:
                    objects[method] = dict.fromkeys(sizes)
                    objects[method][s] = res
                    
                end = time.perf_counter()  
                pr_res, p_dimred = res
                # palantir.plot.plot_palantir_results_f(pr_res, p_dimred)
                # Finer adjustments to the plot can be made using the arguments to this function
                plot_palantir_results_f(pr_res, p_dimred, point_size = 0.3, labels = slabels, origin=origin_new)
                now = datetime.now().strftime('%m-%d-%H%M%S')
                plt.savefig('./plots/'+method+'_'+str(s)+'_pseudotime_'+now)
                plt.close('all')
                
            elif 'paga' in method:
                print('starting paga')
                start = time.perf_counter()
                res = test_paga(sdata, sdimred, labels = labels.iloc[ss], origin=origin_new, ncomps = ncomps)
                dt, oid = res
                if method in objects:
                    objects[method][s] = res
                else:
                    objects[method] = dict.fromkeys(sizes)
                    objects[method][s] = res
                    
                end = time.perf_counter()  
                
                # Finer adjustments to the plot can be made using the arguments to this function
                ax = sc.pl.paga(dt, layout='rt', fontsize=8, fontoutline=1, root=oid, edge_width_scale=0.4, threshold=0.1, show=False)
                handles = create_handles(dt)
                lgd = ax.legend(handles=handles, loc=8, prop={'size': 8}, ncol=len(handles)//5, bbox_to_anchor=(0.5,-0.4))
                now = datetime.now().strftime('%m-%d-%H%M%S')
                plt.savefig('./plots/'+method+'_'+str(s)+'_graph_'+now, bbox_extra_artists=(lgd,), bbox_inches='tight')
                plt.close('all')
                
                # sc.pl.draw_graph(res, color='dpt_pseudotime', layout='rt')
                # plt.savefig('./plots/'+method+'_'+str(s)+'_pseudotime_'+now)
                # plt.close('all')
            
            st += (end-start)
            print("time for iter: {} seconds".format(end-start))
            sys.stdout.flush()

        times[s] = st/repeat
        
    return times, objects
  
def test_palantir(data, dimred, origin = 0, labels = None, n_pca = -1):
    ## preprocessing
    ad = sc.AnnData(data)
    sc.pp.normalize_total(ad) # normalize_per_cell deprecated
    # palantir.preprocess.log_transform(ad)
    # sc.pp.highly_variable_genes(ad, n_top_genes=1500, flavor='cell_ranger')
    if n_pca>0:
        pca_projections, _ = palantir.utils.run_pca(ad, use_hvg=False, n_components=n_pca)
    else:
        pca_projections = pd.DataFrame(ad.X,index=ad.obs_names) #check
        n_pca = ad.X.shape[1]
    # Run diffusion maps
    dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=n_pca)
    # get low dimensional embedding
    ms_data = palantir.utils.determine_multiscale_space(dm_res)
    ## formatting
    dimred.index = ad.obs_names
    
    ## visualize diffusion map
    # palantir.plot.plot_diffusion_components(dimred, dm_res)
    # plt.savefig("diffusion.png")
    print(7)
    
    ## run palantir
    pr_res = palantir.core.run_palantir(ms_data, 
                                        early_cell=ad.obs_names[origin], 
                                        num_waypoints=500)
    return pr_res, dimred
  
def test_via(data, labels = None, ori_id = 0, quality = 'coarse', ncomps=-1):
    
    # setting via parameters
    if ncomps<0:
        ncomps = data.shape[1]-1
    knn=30
    v0_random_seed=4
    root_user = [ori_id] 
    dataset = '' 

    adata = sc.AnnData(data)
    sc.tl.pca(adata, n_comps=ncomps)
    
    v0 = via.VIA(adata.obsm['X_pca'][:, 0:ncomps], labels,
                 jac_std_global=0.15, dist_std_local=1, knn=knn,
                 too_big_factor=0.3, root_user=root_user, dataset=dataset,
                 preserve_disconnected=True, random_seed=v0_random_seed,
                 is_coarse=True,pseudotime_threshold_TS=20,
                 neighboring_terminal_states_threshold=3)

    v0.run_VIA()
    
    if quality=='fine':
        v1 = via.VIA(adata.obsm['X_pca'][:, 0:ncomps], labels, 
                  jac_std_global=0.15, dist_std_local=1, knn=knn, too_big_factor=0.05, 
                  root_user=root_user, x_lazy=0.95, alpha_teleport=0.99, dataset='', 
                  preserve_disconnected=True, super_terminal_clusters=v0.terminal_clusters, 
                  is_coarse=False, random_seed=v0_random_seed, pseudotime_threshold_TS=10, 
                  via_coarse=v0)
        v1.run_VIA()
        return (v0, v1)

    
    return (v0, v0)
  
def test_paga(data, dimred = None, origin = None, labels = None, ncomps = -1):
    # preprocessing
    ad = sc.AnnData(data)
    ad.obs_names_make_unique()
    
    if ncomps>0:
        sc.tl.pca(ad, n_comps=ncomps)

    # Louvain clustering
    sc.pp.neighbors(ad, n_neighbors=30, n_pcs=ncomps)
    sc.tl.louvain(ad, resolution=1) # higher value for more clusters
    ad.obs['clusters'] = ad.obs['louvain']
    
    # run PAGA
    sc.tl.paga(ad, groups='clusters')
    
    # compute pseudotime
    sc.tl.dpt(ad)
    
    # add custom embedding
    ad.obsm['rt'] = dimred.to_numpy()
    
    # match clusters to their most abundant celltype for visualization
    labels.index = ad.obs_names
    corr = match_clusters(ad.obs['clusters'], labels, truncate_labels=-1)
    ad.obs['new_clusters'] = ad.obs['clusters']
    ad.rename_categories('new_clusters', corr)
    
    # find origin cluster
    oname = labels.index[origin]
    oclust = ad.obs['clusters'][oname]
    oid = ad.obs['clusters'].cat.categories.tolist().index(oclust)

    # save data
    ad.write(results_file)
    return ad, oid
  
def match_clusters(clustering, labels, truncate_labels = 5):
    # assuming clusters and labels are pandas df with cell names as index
    cls = list(clustering.cat.categories)
    # correspondances = dict.fromkeys(cls)
    correspondances = []
    for cl in cls:
        wh = clustering[clustering==cl].index
        group = labels.loc[wh]
        freq = group.value_counts()[:2].index.tolist()
        if isinstance(freq[0], tuple):
            mode = freq[0][0] if freq[0][0] != 'ungated' else freq[1][0]
        else:
            mode = freq[0] if freq[0] != 'ungated' else freq[1]
        m = str(mode)[:truncate_labels] if truncate_labels>0 else str(mode)
        correspondances.append(cl+'/'+m)
    return correspondances
  
def create_handles(adata):
    handles = []
    for i, n in enumerate(list(adata.obs.new_clusters.cat.categories)):
        color = adata.uns['clusters_colors'][i]
        label = n.split('/')
        label = ': '.join(label)
        h = Line2D([0], [0], marker='o', color=color, label=label, markerfacecolor=color, markersize=5)
        handles.append(h)
        
    return handles
        
        
        
        
        
        
        
        
        
        
