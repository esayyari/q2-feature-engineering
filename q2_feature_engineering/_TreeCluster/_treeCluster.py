#!/usr/bin/env python3
from q2_types.tree import NewickFormat
from qiime2.plugin import (Int, Float, Str)
import pandas as pd
import biom
import numpy as np
import skbio
from q2_types.tree._transformer import _1, _2
from math import log
from niemads import DisjointSet
from queue import PriorityQueue,Queue
from treeswift import read_tree_newick
from sys import stderr
VERBOSE = False
NUM_THRESH = 1000 # number of thresholds for the threshold-free methods to use

# merge two sorted lists into a sorted list
def merge_two_sorted_lists(x,y):
    out = list(); i = 0; j = 0
    while i < len(x) and j < len(y):
        if x[i] < y[j]:
            out.append(x[i]); i+= 1
        else:
            out.append(y[j]); j += 1
    while i < len(x):
        out.append(x[i]); i += 1
    while j < len(y):
        out.append(y[j]); j += 1
    return out

# merge multiple sorted lists into a sorted list
def merge_multi_sorted_lists(lists):
    pq = PriorityQueue()
    for l in range(len(lists)):
        if len(lists[l]) != 0:
            pq.put((lists[l][0],l))
    inds = [1 for _ in range(len(lists))]
    out = list()
    while not pq.empty():
        d,l = pq.get(); out.append(d)
        if inds[l] < len(lists[l]):
            pq.put((lists[l][inds[l]],l)); l += 1
    return out

# get the median of a sorted list
def median(x):
    if len(x) % 2 != 0:
        return x[int(len(x)/2)]
    else:
        return (x[int(len(x)/2)]+x[int(len(x)/2)-1])/2

# get the average of a list
def avg(x):
    return float(sum(x))/len(x)

# convert p-distance to Jukes-Cantor distance
def p_to_jc(d,seq_type):
    b = {'dna':3./4., 'protein':19./20.}[seq_type]
    return -1*b*log(1-(d/b))

# cut out the current node's subtree (by setting all nodes' DELETED to True) and return list of leaves
def cut(node):
    cluster = list()
    descendants = Queue(); descendants.put(node)
    while not descendants.empty():
        descendant = descendants.get()
        if descendant.DELETED:
            continue
        descendant.DELETED = True
        descendant.left_dist = 0; descendant.right_dist = 0; descendant.edge_length = 0
        if descendant.is_leaf():
            cluster.append(str(descendant))
        else:
            for c in descendant.children:
                descendants.put(c)
    return cluster

# initialize properties of input tree and return set containing taxa of leaves
def prep(tree,support):
    tree.resolve_polytomies(); tree.suppress_unifurcations()
    leaves = set()
    for node in tree.traverse_postorder():
        if node.edge_length is None:
            node.edge_length = 0
        node.DELETED = False
        if node.is_leaf():
            leaves.add(str(node))
        else:
            try:
                node.confidence = float(str(node))
            except:
                node.confidence = 100. # give edges without support values support 100
            if node.confidence < support: # don't allow low-support edges
                node.edge_length = float('inf')
    return leaves

# return a sorted list of all unique pairwise leaf distances <= a given threshold
def pairwise_dists_below_thresh(tree,threshold):
    pairwise_dists = set()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.leaf_dists = {0}; node.min_leaf_dist = 0
        else:
            children = list(node.children)
            for i in range(len(children)-1):
                c1 = children[i]
                for j in range(i+1,len(children)):
                    c2 = children[j]
                    for d1 in c1.leaf_dists:
                        for d2 in c2.leaf_dists:
                            pd = d1 + c1.edge_length + d2 + c2.edge_length
                            if pd <= threshold:
                                pairwise_dists.add(pd)
            node.leaf_dists = set(); node.min_leaf_dist = float('inf')
            for c in children:
                if c.min_leaf_dist + c.edge_length > threshold:
                    continue
                for d in c.leaf_dists:
                    nd = d+c.edge_length
                    if nd < threshold:
                        node.leaf_dists.add(nd)
                    if nd < node.min_leaf_dist:
                        node.min_leaf_dist = nd
    return sorted(pairwise_dists)

# split leaves into minimum number of clusters such that the maximum leaf pairwise distance is below some threshold
def min_clusters_threshold_max(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        if node.is_leaf():
            node.left_dist = 0; node.right_dist = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node); continue
            if children[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(children[0].left_dist,children[0].right_dist) + children[0].edge_length
            if children[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(children[1].left_dist,children[1].right_dist) + children[1].edge_length

            # if my kids are screwing things up, cut out the longer one
            if node.left_dist + node.right_dist > threshold:
                if node.left_dist > node.right_dist:
                    cluster = cut(children[0])
                    node.left_dist = 0
                else:
                    cluster = cut(children[1])
                    node.right_dist = 0

                # add cluster
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# median leaf pairwise distance cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_med_clade(tree,threshold,support):
    leaves = prep(tree,support)

    # bottom-up traversal to compute median pairwise distances
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.med_pair_dist = 0
            node.leaf_dists = [0]
            node.pair_dists = list()
        else:
            children = list(node.children)
            l_leaf_dists = [d + children[0].edge_length for d in children[0].leaf_dists]
            r_leaf_dists = [d + children[1].edge_length for d in children[1].leaf_dists]
            node.leaf_dists = merge_two_sorted_lists(l_leaf_dists,r_leaf_dists)
            if len(l_leaf_dists) < len(r_leaf_dists):
                across_leaf_dists = [[l+r for r in r_leaf_dists] for l in l_leaf_dists]
            else:
                across_leaf_dists = [[l+r for l in l_leaf_dists] for r in r_leaf_dists]
            node.pair_dists = merge_multi_sorted_lists([children[0].pair_dists,children[1].pair_dists] + across_leaf_dists)
            if node.pair_dists[-1] == float('inf'):
                node.med_pair_dist = float('inf')
            else:
                node.med_pair_dist = median(node.pair_dists)
            for c in (children[0],children[1]):
                del c.leaf_dists; del c.pair_dists

    # perform clustering
    q = Queue(); q.put(tree.root); roots = list()
    while not q.empty():
        node = q.get()
        if node.med_pair_dist <= threshold:
            roots.append(node)
        else:
            for c in node.children:
                q.put(c)

    # if verbose, print the clades defined by each cluster
    if VERBOSE:
        for root in roots:
            print("%s;" % root.newick(), file=stderr)
    return [[str(l) for l in root.traverse_leaves()] for root in roots]

# average leaf pairwise distance cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_avg_clade(tree,threshold,support):
    leaves = prep(tree,support)

    # bottom-up traversal to compute average pairwise distances
    for node in tree.traverse_postorder():
        node.total_pair_dist = 0; node.total_leaf_dist = 0
        if node.is_leaf():
            node.num_leaves = 1
            node.avg_pair_dist = 0
        else:
            children = list(node.children)
            node.num_leaves = sum(c.num_leaves for c in children)
            node.total_pair_dist = children[0].total_pair_dist + children[1].total_pair_dist + (children[0].total_leaf_dist*children[1].num_leaves + children[1].total_leaf_dist*children[0].num_leaves)
            node.total_leaf_dist = (children[0].total_leaf_dist + children[0].edge_length*children[0].num_leaves) + (children[1].total_leaf_dist + children[1].edge_length*children[1].num_leaves)
            node.avg_pair_dist = node.total_pair_dist/((node.num_leaves*(node.num_leaves-1))/2)

    # perform clustering
    q = Queue(); q.put(tree.root); roots = list()
    while not q.empty():
        node = q.get()
        if node.avg_pair_dist <= threshold:
            roots.append(node)
        else:
            for c in node.children:
                q.put(c)

    # if verbose, print the clades defined by each cluster
    if VERBOSE:
        for root in roots:
            print("%s;" % root.newick(), file=stderr)
    return [[str(l) for l in root.traverse_leaves()] for root in roots]

# total branch length cannot exceed threshold, and clusters must define clades
def min_clusters_threshold_sum_bl_clade(tree,threshold,support):
    leaves = prep(tree,support)

    # compute branch length sums
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.total_bl = 0
        else:
            node.total_bl = sum(c.total_bl + c.edge_length for c in node.children)

    # perform clustering
    q = Queue(); q.put(tree.root); roots = list()
    while not q.empty():
        node = q.get()
        if node.total_bl <= threshold:
            roots.append(node)
        else:
            for c in node.children:
                q.put(c)

    # if verbose, print the clades defined by each cluster
    if VERBOSE:
        for root in roots:
            print("%s;" % root.newick(), file=stderr)
    return [[str(l) for l in root.traverse_leaves()] for root in roots]

# total branch length cannot exceed threshold
def min_clusters_threshold_sum_bl(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.left_total = 0; node.right_total = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node); continue
            if children[0].DELETED:
                node.left_total = 0
            else:
                node.left_total = children[0].left_total + children[0].right_total + children[0].edge_length
            if children[1].DELETED:
                node.right_total = 0
            else:
                node.right_total = children[1].left_total + children[1].right_total + children[1].edge_length
            if node.left_total + node.right_total > threshold:
                if node.left_total > node.right_total:
                    cluster = cut(children[0])
                    node.left_total = 0
                else:
                    cluster = cut(children[1])
                    node.right_total = 0
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# single-linkage clustering using Metin's cut algorithm
def single_linkage_cut(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    # find closest leaf below (dist,leaf)
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.min_below = (0,node.label)
        else:
            node.min_below = min((c.min_below[0]+c.edge_length,c.min_below[1]) for c in node.children)

    # find closest leaf above (dist,leaf)
    for node in tree.traverse_preorder():
        node.min_above = (float('inf'),None)
        if node.is_root():
            continue
        # min distance through sibling
        for c in node.parent.children:
            if c != node:
                dist = node.edge_length + c.edge_length + c.min_below[0]
                if dist < node.min_above[0]:
                    node.min_above = (dist,c.min_below[1])
        # min distance through grandparent
        if not c.parent.is_root():
            dist = node.edge_length + node.parent.min_above[0]
            if dist < node.min_above[0]:
                node.min_above = (dist,node.parent.min_above[1])

    # find clusters
    for node in tree.traverse_postorder(leaves=False):
        # assume binary tree here (prep function guarantees this)
        l_child,r_child = node.children
        l_dist = l_child.min_below[0] + l_child.edge_length
        r_dist = r_child.min_below[0] + r_child.edge_length
        a_dist = node.min_above[0]
        bad = [0,0,0] # left, right, up
        if l_dist + r_dist > threshold:
            bad[0] += 1; bad[1] += 1
        if l_dist + a_dist > threshold:
            bad[0] += 1; bad[2] += 1
        if r_dist + a_dist > threshold:
            bad[1] += 1; bad[2] += 1
        # cut either (or both) children
        for i in [0,1]:
            if bad[i] == 2:
                cluster = cut(node.children[i])
                if len(cluster) != 0:
                    clusters.append(cluster)
                    for leaf in cluster:
                        leaves.remove(leaf)
        # cut above (equals cutting me)
        if bad[2] == 2: # if cutting above, just cut me
            cluster = cut(node)
            if len(cluster) != 0:
                clusters.append(cluster)
                for leaf in cluster:
                    leaves.remove(leaf)
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# single-linkage clustering using Niema's union algorithm
def single_linkage_union(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()

    # find closest leaf below (dist,leaf)
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.min_below = (0,node.label)
        else:
            node.min_below = min((c.min_below[0]+c.edge_length,c.min_below[1]) for c in node.children)

    # find closest leaf above (dist,leaf)
    for node in tree.traverse_preorder():
        node.min_above = (float('inf'),None)
        if node.is_root():
            continue
        # min distance through sibling
        for c in node.parent.children:
            if c != node:
                dist = node.edge_length + c.edge_length + c.min_below[0]
                if dist < node.min_above[0]:
                    node.min_above = (dist,c.min_below[1])
        # min distance through grandparent
        if not c.parent.is_root():
            dist = node.edge_length + node.parent.min_above[0]
            if dist < node.min_above[0]:
                node.min_above = (dist,node.parent.min_above[1])

    # set up Disjoint Set
    ds = DisjointSet(leaves)
    for node in tree.traverse_preorder(leaves=False):
        # children to min above
        for c in node.children:
            if c.min_below[0] + c.edge_length + node.min_above[0] <= threshold:
                ds.union(c.min_below[1], node.min_above[1])
        for i in range(len(node.children)-1):
            c1 = node.children[i]
            for j in range(i+1, len(node.children)):
                c2 = node.children[j]
                if c1.min_below[0] + c1.edge_length + c2.min_below[0] + c2.edge_length <= threshold:
                    ds.union(c1.min_below[1], c2.min_below[1])
    return [list(s) for s in ds.sets()]

# min_clusters_threshold_max, but all clusters must define a clade
def min_clusters_threshold_max_clade(tree,threshold,support):
    leaves = prep(tree,support)

    # compute leaf distances and max pairwise distances
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.leaf_dist = 0
            node.max_pair_dist = 0
        else:
            node.leaf_dist = max(c.leaf_dist + c.edge_length for c in node.children)
            curr_pair_dist = sum(c.leaf_dist + c.edge_length for c in node.children) # bifurcating because of prep
            node.max_pair_dist = max([c.max_pair_dist for c in node.children] + [curr_pair_dist])

    # perform clustering
    q = Queue(); q.put(tree.root); roots = list()
    while not q.empty():
        node = q.get()
        if node.max_pair_dist <= threshold:
            roots.append(node)
        else:
            for c in node.children:
                q.put(c)

    # if verbose, print the clades defined by each cluster
    if VERBOSE:
        for root in roots:
            print("%s;" % root.newick(), file=stderr)
    return [[str(l) for l in root.traverse_leaves()] for root in roots]

# pick the threshold between 0 and "threshold" that maximizes number of (non-singleton) clusters
def argmax_clusters(method,tree,threshold,support):
    from copy import deepcopy
    assert threshold > 0, "Threshold must be positive"
    #thresholds = pairwise_dists_below_thresh(deepcopy(tree),threshold)
    thresholds = [i*threshold/NUM_THRESH for i in range(NUM_THRESH+1)]
    best = None; best_num = -1; best_t = -1
    for i,t in enumerate(thresholds):
        if VERBOSE:
            print("%s%%"%str(i*100/len(thresholds)).rstrip('0'),end='\r',file=stderr)
        clusters = method(deepcopy(tree),t,support)
        num_non_singleton = len([c for c in clusters if len(c) > 1])
        if num_non_singleton > best_num:
            best = clusters; best_num = num_non_singleton; best_t = t
    print("\nBest Threshold: %f"%best_t,file=stderr)
    return best

# cut all branches longer than the threshold
def length(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_postorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # if i'm screwing things up, cut me
        if node.edge_length is not None and node.edge_length > threshold:
            cluster = cut(node)
            if len(cluster) != 0:
                clusters.append(cluster)
                for leaf in cluster:
                    leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# same as length, and clusters must define a clade
def length_clade(tree,threshold,support):
    leaves = prep(tree,support)

    # compute max branch length in clades
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.max_bl = 0
        else:
            node.max_bl = max([c.max_bl for c in node.children] + [c.edge_length for c in node.children])

    # perform clustering
    q = Queue(); q.put(tree.root); roots = list()
    while not q.empty():
        node = q.get()
        if node.max_bl <= threshold:
            roots.append(node)
        else:
            for c in node.children:
                q.put(c)

    # if verbose, print the clades defined by each cluster
    if VERBOSE:
        for root in roots:
            print("%s;" % root.newick(), file=stderr)
    return [[str(l) for l in root.traverse_leaves()] for root in roots]

# cut tree at threshold distance from root (clusters will be clades by definition) (ignores support threshold if branch is below cutting point)
def root_dist(tree,threshold,support):
    leaves = prep(tree,support)
    clusters = list()
    for node in tree.traverse_preorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue
        if node.is_root():
            node.root_dist = 0
        else:
            node.root_dist = node.parent.root_dist + node.edge_length
        if node.root_dist > threshold:
            cluster = cut(node)
            if len(cluster) != 0:
                clusters.append(cluster)
                for leaf in cluster:
                    leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters

# cut tree at threshold distance from the leaves (if tree not ultrametric, max = distance from furthest leaf from root, min = distance from closest leaf to root, avg = average of all leaves)
def leaf_dist(tree,threshold,support,mode):
    modes = {'max':max,'min':min,'avg':avg}
    assert mode in modes, "Invalid mode. Must be one of: %s" % ', '.join(sorted(modes.keys()))
    dist_from_root = modes[mode](d for u,d in tree.distances_from_root(internal=False)) - threshold
    return root_dist(tree,dist_from_root,support)
def leaf_dist_max(tree,threshold,support):
    return leaf_dist(tree,threshold,support,'max')
def leaf_dist_min(tree,threshold,support):
    return leaf_dist(tree,threshold,support,'min')
def leaf_dist_avg(tree,threshold,support):
    return leaf_dist(tree,threshold,support,'avg')

METHODS = {
    'max': min_clusters_threshold_max,
    'max_clade': min_clusters_threshold_max_clade,
    'sum_branch': min_clusters_threshold_sum_bl,
    'sum_branch_clade': min_clusters_threshold_sum_bl_clade,
    'avg_clade': min_clusters_threshold_avg_clade,
    'med_clade': min_clusters_threshold_med_clade,
    'single_linkage': single_linkage_cut,
    'single_linkage_cut': single_linkage_cut,
    'single_linkage_union': single_linkage_union,
    'length': length,
    'length_clade': length_clade,
    'root_dist': root_dist,
    'leaf_dist_max': leaf_dist_max,
    'leaf_dist_min': leaf_dist_min,
    'leaf_dist_avg': leaf_dist_avg
}
THRESHOLDFREE = {'argmax_clusters':argmax_clusters}


def _clean_phylogeny(phylogeny_fp, table):
    phylogeny = _2(phylogeny_fp)
    obs = list(table.ids('observation'))
    phylogeny = phylogeny.shear(obs)
    phylogeny.bifurcate()
    for node in phylogeny.traverse():
        if node.length is None:
            node.length = 0
    tree = phylogeny.copy()
    tree_pruned_out = _1(tree)
    return tree_pruned_out


def _read_phylogeny(phylogeny_fp):

    with open(str(phylogeny_fp)) as input_file:
        tree_str = input_file.readline()
        if isinstance(tree_str, bytes):
            tree_str_decoded = tree_str.decode().strip()
        else:
            tree_str_decoded = tree_str.strip()
        tree = read_tree_newick(tree_str_decoded)

    return tree


def initialize_mapping_matrix(mapping, table):
    m = max(mapping.Taxon)
    n = len(table.ids('observation'))
    k = max(mapping.Taxon)
    if m < 0 or n < 0:
        return []
    map_matrix = np.zeros((m, n))
    return map_matrix


def pivot_mapping(mapping_df):
    clst_dict = {}
    for clst in mapping_df['Taxon'].unique():
        clst_dict[clst] = mapping_df.loc[mapping_df.Taxon == clst].index.values
    return clst_dict


def make_mapping_matrix(mapping, map_matrix):
    i = 0
    j = 0
    obs = [""] * len(map_matrix)
    singletones = []
    for idx, row in mapping.iterrows():
        if row['Taxon'] == -1:
            singletones.append(i)
            j += 1
            obs.append(idx)
        else:
            map_matrix[row['Taxon']-1][i] = 1
            obs[row['Taxon']-1] = 'cluster' + str(row['Taxon'])
        i += 1
    singletones = np.asarray(singletones)
    pivot_mapping_dct = pivot_mapping(mapping)
    return map_matrix, singletones, obs, pivot_mapping_dct


def check_mapping(final_data, table, new_data, mapping, map_matrix):
    flag = False
    check = []

    if not final_data.shape[0] == new_data.shape[0] + (mapping.Taxon == -1).sum():
        print("The number of features in the mapped and original tables are inconsistent!")
        print(final_data.shape[0], new_data.shape[0] + (mapping.Taxon == -1).sum())
        flag = True
    if not final_data.sum() == table.matrix_data.sum():
        print("The sum of all faetures in the original table and mapped one are not the same!")
        print(final_data.sum(), table.matrix_data.sum())
        flag = True
    for i in range(new_data.shape[0]):
        check.append(set(mapping.iloc[np.where(map_matrix[i,:])].Taxon.values) == set([i+1]))
    if np.sum(check) != len(check):
        print("There is an inconsistency in the mapped features")
        flag = True
    if not flag:
        print("The feature mapping finished successfully! ")
    return flag


def map_features(mapping, table):
    mapping = mapping.reset_index().set_index('Feature ID')
    mapping = mapping.loc[table.ids('observation')].reindex(table.ids('observation'))
    map_matrix = initialize_mapping_matrix(mapping, table)
    map_matrix, singletones, obs, pivot_mapping_dct = make_mapping_matrix(mapping, map_matrix)
    if len(map_matrix) > 0:
        new_data = map_matrix * table.matrix_data
        final_data = np.row_stack([new_data, np.asarray(table.matrix_data[singletones, :].todense())])

    else:
        final_data = np.asarray(table.matrix_data[singletones, :].todense())

    mapped_biom_table = biom.Table(data=final_data, observation_ids=obs, sample_ids=table.ids())
    if len(map_matrix) > 0 and check_mapping(final_data, table, new_data, mapping, map_matrix):
        raise ValueError("something went wrong!")
    return mapped_biom_table, pivot_mapping_dct


def _convert_result_to_data_frame(clusters):
    cluster_num = 1
    frag_name = []
    clust_name = []
    for cluster in clusters:
        if len(cluster) == 1:
            frag_name.append(list(cluster)[0])
            clust_name.append(-1)
        else:
            for sp in cluster:
                frag_name.append(sp)
                clust_name.append(cluster_num)
            cluster_num += 1

    cluster_dct = {'Feature ID': frag_name, 'Taxon': clust_name}
    cluster_df = pd.DataFrame(cluster_dct)
    return cluster_df


def average_branch_lengths(phylogeny_fp):
    tree = _2(phylogeny_fp)
    j = 0
    for nd in tree.postorder():
        if nd.is_tip():
            nd.num = 1
            nd.avg = 0.0
            nd.sum_dist = 0.0
            nd.num_pairs = 1
            nd.dummy_name = nd.name
        else:
            child_nodes = nd.children
            nd.num_pairs = child_nodes[0].num * child_nodes[1].num
            total_dist = child_nodes[1].sum_dist * child_nodes[0].num + child_nodes[0].sum_dist * child_nodes[1].num + (
                        child_nodes[1].length + child_nodes[0].length) * nd.num_pairs
            nd.avg = total_dist / nd.num_pairs / 2.0

            nd.sum_dist = child_nodes[0].sum_dist + child_nodes[1].sum_dist + child_nodes[0].length * child_nodes[
                0].num + child_nodes[1].length * child_nodes[1].num

            nd.num = child_nodes[0].num + child_nodes[1].num
            nd.dummy_name = "dummy" + str(j)
            j += 1
    return tree


def generate_new_phylogeny(phylogeny_fp, pivot_mapping_dct):
    tree = average_branch_lengths(phylogeny_fp)
    observed_internal_nodes = set()
    for idx in pivot_mapping_dct.keys():
        if idx == -1:
            continue
        else:
            node_lca = tree.lca(pivot_mapping_dct[idx])
            if len(list(node_lca.tips())) != len(pivot_mapping_dct[idx]):
                raise ValueError("The LCA has different tips than expected", len(list(node_lca.tips())),
                                 len(pivot_mapping_dct[idx]))
            tmp_node = skbio.TreeNode('cluster' + str(idx))
            tmp_node.length = node_lca.avg
            parent_node = node_lca.parent
            if node_lca.dummy_name in observed_internal_nodes:
                raise ValueError("node_lca was observed before!", node_lca.dummy_name, node_lca, idx)
            else:
                observed_internal_nodes.add(node_lca.dummy_name)
            parent_node.append(tmp_node)
            parent_node.remove(node_lca)
    return tree


def tree_cluster(phylogeny_fp: NewickFormat, table: biom.Table, method: Str = 'max_clade', threshold: Float = 0.0,
                 threshold_free: Str = None) -> (skbio.TreeNode, pd.DataFrame, biom.Table):
    phylogeny_cleaned_fp = _clean_phylogeny(phylogeny_fp, table)
    print("The numbe of features in the original table is", len(table.ids('observation')))

    print("The phylogeny loaded and randomly resolved!")
    tree = _read_phylogeny(str(phylogeny_cleaned_fp))
    if threshold_free is None:
        clusters = METHODS[method.lower()](tree, threshold, float('-inf'))
    else:
        clusters = THRESHOLDFREE[threshold_free](METHODS[method.lower()], tree, threshold, float('-inf'))

    cluster_df = _convert_result_to_data_frame(clusters)
    print("TreeCluster finished successfully! The number of clusters using the threshold",
          threshold, "is", len(clusters))

    final_table_filtered, pivot_mapping_dct = map_features(cluster_df, table)
    print("The clustered table is generated!")
    new_tree = generate_new_phylogeny(phylogeny_cleaned_fp, pivot_mapping_dct)
    print("The clustered phylogeny is generated!")
    cluster_df = cluster_df.set_index('Feature ID')
    return new_tree, cluster_df, final_table_filtered
