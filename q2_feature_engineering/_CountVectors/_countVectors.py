from skbio.diversity._util import _vectorize_counts_and_tree
from skbio.diversity.beta._unifrac import _setup_multiple_unweighted_unifrac, \
    _setup_multiple_weighted_unifrac, \
    _get_tip_indices
import biom
import numpy as np
import skbio
from time import time
import tempfile


def _map_observations(table: biom.Table) -> biom.Table:
    obs_dict = {}
    for taxa in table.ids('observation'):
        obs_dict[taxa] = taxa.replace('_', ' ')
    table = table.update_ids(id_map=obs_dict, axis='observation', inplace=False)
    return table


def rename_otus(tree_index):
    new_otus = tree_index['name']
    j = 0
    for i, name in enumerate(new_otus):
        if not name:
            new_otus[i] = "None_" + str(j)
        else:
            new_otus[i] += "_" + str(j)
        j += 1
    return new_otus


def prune_features_from_phylogeny(table: biom.Table,
                                  phylogeny: skbio.TreeNode) -> skbio.TreeNode:
    print('Will prune the phylogeny')
    tree = phylogeny
    obs = table.ids('observation')
    tip_names_set = set([x.name for x in tree.tips()])
    to_delete_names = tip_names_set - set(obs)
    to_delete_set = to_delete_names

    if len(set(obs) - tip_names_set) > 0:
        raise ValueError(
            "There are",  len(set(obs) - tip_names_set),
            "features in the feature table not present in the phylogeny! "
            "Please check your tree"
        )
    else:
        print("All", len(obs), "features present in the feature table are also "
                               "in the phylogeny.")

    if len(to_delete_set) > 0:
        t0 = time()
        print("The set of features in the phylogeny and the table "
              "are not the same.",
              len(to_delete_set),
              "features will be pruned from the tree.")
        tree_pruned = tree.shear(set(obs))
        print("It takes", time()-t0, "seconds to prune the phylogeny")
        to_delete_set = set([x.name for x in tree_pruned.tips()]) - set(obs)
        to_delete_rev_set = set(obs) - set([x.name for x in tree_pruned.tips()])
        if len(to_delete_set) > 0 or len(to_delete_rev_set):
            raise ValueError(
                "Pruning the phylogeny failed! There are", len(to_delete_set),
                "features in the phylogeny not "
                "present in the feature table, and",
                len(to_delete_rev_set),
                "features in the feature table not available in the "
                "phylogeny! Both should be 0"
            )
        else:
            print("The phylogeny was pruned successfully!")
    else:
        print("The set of features in the phylogeny and "
              "the table are the same. "
              "No feature will be pruned from the tree.")
        tree_pruned = tree
    return tree_pruned


def count_vectors(table: biom.Table, phylogeny: skbio.TreeNode,
                        method: str = 'weighted_unifrac') -> biom.Table:
    table = _map_observations(table)
    pruned_phylo = prune_features_from_phylogeny(table, phylogeny)
    pruned_phylo = rename_nodes(pruned_phylo)
    table = table.sort(axis='observation')
    otu_ids = np.asarray(table.ids('observation'))
    counts = np.asarray(table.matrix_data.todense().transpose())
    features, tree_index = _run_unifrac(counts, otu_ids, pruned_phylo, method)
    return biom.Table(data=features.transpose(),
                      observation_ids=rename_otus(tree_index),
                      sample_ids=table.ids())


def rename_nodes(phylogeny):
    for nd in phylogeny.postorder():
        if nd.is_tip():
            nd.list_nodes = [nd.name]
        else:
            nd.list_nodes = []
            for ch in nd.children:
                nd.list_nodes += ch.list_nodes
            nd.name = "_".join(sorted(nd.list_nodes))
    return phylogeny


def _run_unifrac(counts, otu_ids, pruned_phylo, method):
    _, tree_index, branch_lengths = \
        _vectorize_counts_and_tree(counts[0, :], otu_ids, pruned_phylo)
    if method == 'weighted_unifrac':
        features = _weighted_unifrac_features(counts,
                                              otu_ids, pruned_phylo,
                                              tree_index, branch_lengths)
    elif method == 'unweighted_unifrac':
        features = _unweighted_unifrac_features(counts,
                                                branch_lengths,
                                                otu_ids,
                                                pruned_phylo)
    else:
        raise ValueError('Method not implemented. Options are '
                         'weighted_unifrac or unweighted_unifrac, but',
                         method, 'is given')
    return features, tree_index


def _get_total_counts(counts, tip_indices):
    counts_sum_list = []
    for i in range(len(counts)):
        counts_sum_list.append(np.take(counts[i], tip_indices).sum())
    return np.asarray(counts_sum_list).reshape((len(counts_sum_list), 1))


def _weighted_unifrac_features(counts, otu_ids, pruned_phylo,
                               tree_index, branch_lengths):
    counts = np.asarray(counts)
    _, ordered_counts_by_node = \
        _setup_multiple_weighted_unifrac(counts, otu_ids, pruned_phylo,
                                         True, True)
    tip_indices = _get_tip_indices(tree_index)
    counts_sum_arr = _get_total_counts(ordered_counts_by_node, tip_indices)
    ordered_counts_by_node = ordered_counts_by_node / counts_sum_arr
    features = np.multiply(branch_lengths, ordered_counts_by_node)
    return features


def _unweighted_unifrac_features(counts, branch_lengths, otu_ids, pruned_phylo):
    _, ordered_counts_by_node = \
        _setup_multiple_unweighted_unifrac(counts, otu_ids, pruned_phylo, True)
    _normalized_counts = np.logical_xor(0, ordered_counts_by_node)
    features = np.multiply(_normalized_counts, branch_lengths)
    return features



