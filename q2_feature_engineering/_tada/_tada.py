# ----------------------------------------------------------------------------
# Copyright (c) 2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from skbio.tree import TreeNode
import skbio
from qiime2.plugin import (Str, Int, Bool, Float, Metadata)
from q2_types.tree import NewickFormat
from qiime2 import NumericMetadataColumn
from .TADA_utils import *
import pandas as pd
import biom
import tempfile
from .datasets import make_data_frame
import shutil
import qiime2
import numpy as np
from q2_types.tree._transformer import _1
from qiime2 import Artifact


def _map_observations(table: biom.Table) -> biom.Table:
    obs_dict = {}
    for taxa in table.ids('observation'):
        obs_dict[taxa] = taxa.replace('_', ' ')
    table = table.update_ids(id_map=obs_dict, axis='observation', inplace=False)
    return table


def _sort_metada(targets_metadata, biom_table):
    targets = targets_metadata.to_dataframe()
    # filter features and targest so samples match
    index = set(targets.index)
    index = [ix for ix in biom_table.ids('sample') if ix in index]
    targets = targets.loc[index]
    feature_data = biom_table.filter(index, inplace=False)

    samples_sorted = np.sort(feature_data.ids('sample'))
    feature_data = feature_data.sort_order(samples_sorted)
    targets = targets.reindex(samples_sorted)
    return targets, feature_data


def _read_inputs(biom_table: biom.Table, phylogeny_fp: NewickFormat,
                 meta_data: NumericMetadataColumn = None):
    if meta_data:
        generate_strategy = "balancing"
        meta, biom_table = _sort_metada(meta_data, biom_table)
        y = meta.iloc[:, 0]
        samples = meta.index
    else:
        generate_strategy = "augmentation"
        y = pd.Series(data=np.ones((len(biom_table.ids('sample')),)),
                      index=biom_table.ids('sample'))
        samples = biom_table.ids('sample')

    _table_tmp = biom_table.sort_order(axis='sample', order=samples)
    _table = _map_observations(_table_tmp)
    pruned_phylogeny_fp = _prune_features_from_phylogeny(_table, phylogeny_fp)
    _tree = dendropy.Tree.get(path=str(pruned_phylogeny_fp),
                              preserve_underscores=False,
                              schema="newick", rooting='default-rooted')

    if sum(samples != _table.ids('sample')) > 0:
        raise ValueError("The samples IDs in meta data and biom table are "
                         "not the same! The difference is:",
                         set(samples) - set(_table.ids('sample')),
                         "Please double check.")

    return _table, y, _tree, generate_strategy, pruned_phylogeny_fp


def tada(phylogeny: NewickFormat, table: biom.Table,
         meta_data: NumericMetadataColumn = None,
         seed_num: Int = 0, xgen: Int = 0,
         n_beta: Int = 1, n_binom: Int = 5, var_method: Str = 'br_penalized',
         stat_method: Str = 'binom', prior_weight: Float = 0,
         coef: Float = 200, exponent: Float = 0.5,
         pseudo_branch_length: Float = 1e-6, pseudo_cnt: Float = 5,
         normalized: Bool = False, output_log_fp: Str = None,
         original_table: Str = None, augmented_table: Str = None,
         concatenate_meta: Metadata = None,
         sampling_strategy: Str = None) -> (NewickFormat, biom.Table):
    _table, y, _phylogeny, generate_strategy, pruned_phylogeny = \
        _read_inputs(biom_table=table, phylogeny_fp=phylogeny,
                     meta_data=meta_data)
    if generate_strategy is 'balancing' and (concatenate_meta is None):
        raise ValueError(
            "Expected a path to write out the generated and original labels and metadata!"
        )
    tmp = tempfile.mkdtemp()
    try:
        sG = SampleGenerator(seed_num=seed_num, logger=None,
                             generate_strategy=generate_strategy, tmp_dir=tmp,
                             xgen=xgen, n_beta=n_beta, n_binom=n_binom,
                             var_method=var_method, stat_method=stat_method,
                             prior_weight=prior_weight,
                             coef=coef, exponent=exponent,
                             pseudo=pseudo_branch_length,
                             pseudo_cnt=pseudo_cnt, normalized=normalized)

        orig_biom, orig_labels, augm_biom, augm_labels = \
            sG.fit_transform(table=_table, y=y, tree=_phylogeny,
                             sampling_strategy=sampling_strategy)
        if np.sum(orig_biom.matrix_data - table.matrix_data) > 1e-20:
            raise ValueError(
                "The original biom table doesn't match the "
                "output of generator function! Please double check")
        if generate_strategy is 'balancing':
            orig_pd, augm_pd = make_data_frame(orig_biom, augm_biom, orig_labels, augm_labels)

            if concatenate_meta and not os.path.exists(os.path.dirname(str(concatenate_meta))):
                os.mkdir(os.path.dirname(str(concatenate_meta)))

            concat_pd = pd.concat([orig_pd, augm_pd])
            concat_meta = qiime2.Metadata(concat_pd)
            concat_meta.save(concatenate_meta)

        if output_log_fp is not None:
            if not os.path.exists(os.path.dirname(output_log_fp)):
                os.mkdir(os.path.dirname(output_log_fp))
            shutil.copyfile(sG.log_fp, output_log_fp)
        if np.sum(orig_biom.ids('observation') == augm_biom.ids('observation'))\
                != len(orig_biom.ids('observation')):
            raise ValueError(
                "The order of features in original and augmented data "
                "is different. Please make sure that your phylogeny doesn't "
                "have extra features"
            )
        if original_table and \
                not os.path.exists(os.path.dirname(str(original_table))):
            os.mkdir(os.path.dirname(str(original_table)))
        elif augmented_table and \
                not os.path.exists(os.path.dirname(str(original_table))):
            os.mkdir(os.path.dirname(str(augmented_table)))
        if augmented_table is not None:
            augm_qza = \
                Artifact.import_data("FeatureTable[Frequency]", augm_biom)
            augm_qza.save(augmented_table)
        if original_table is not None:
            orig_qza = \
                Artifact.import_data("FeatureTable[Frequency]", orig_biom)
            orig_qza.save(original_table)

        concat_biom = orig_biom
        concat_biom = concat_biom.merge(augm_biom)
    finally:
        print("Something went wrong")

    return pruned_phylogeny, concat_biom


def prune_features_from_phylogeny(table: biom.Table,
                                  phylogeny: skbio.TreeNode) -> skbio.TreeNode:
    print('Will prune the phylogeny')
    # tree = TreeNode.read(str(phylogeny))
    tree = phylogeny
    table = _map_observations(table)
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
        print("All", len(obs), "features present in the feature table are "
                               "also in the phylogeny.")

    if len(to_delete_set) > 0:
        t0 = time()
        print("The set of features in the phylogeny and the table are "
              "not the same.", len(to_delete_set), "features will be "
                                                   "pruned from the tree.")
        tree_pruned = tree.shear(set(obs))
        print("It takes", time()-t0, "seconds to prune the phylogeny")
        to_delete_set = set([x.name for x in tree_pruned.tips()]) - set(obs)
        to_delete_rev_set = set(obs) - set([x.name for x in tree_pruned.tips()])
        if len(to_delete_set) > 0 or len(to_delete_rev_set):
            raise ValueError(
                "Pruning the phylogeny failed! There are",
                len(to_delete_set), "features in the phylogeny not "
                                    "present in the feature table, and",
                len(to_delete_rev_set), "features in the feature table "
                                        "not available in the phylogeny! "
                                        "Both should be 0"
            )
        else:
            print("The phylogeny was pruned successfully!")
    else:
        print("The set of features in the phylogeny and the "
              "table are the same. No feature will be pruned from the tree.")
        tree_pruned = tree
    return tree_pruned


def _prune_features_from_phylogeny(table: biom.Table,
                                   phylogeny_fp: NewickFormat) -> NewickFormat:
    print('Will prune the phylogeny')
    tree = TreeNode.read(str(phylogeny_fp))
    obs = table.ids('observation')
    tip_names_set = set([x.name for x in tree.tips()])
    to_delete_names = tip_names_set - set(obs)
    to_delete_set = to_delete_names

    if len(set(obs) - tip_names_set) > 0:
        raise ValueError(
            "There are",  len(set(obs) - tip_names_set),
            "features in the feature table not present "
            "in the phylogeny! Please check your tree"
        )
    else:
        print("All", len(obs), "features present in the "
                               "feature table are also in the phylogeny.")

    if len(to_delete_set) > 0:
        t0 = time()
        print("The set of features in the phylogeny and the table "
              "are not the same.", len(to_delete_set),
              "features will be pruned from the tree.")
        tree_pruned = tree.shear(set(obs))
        print("It takes", time()-t0, "seconds to prune the phylogeny")
        to_delete_set = set([x.name for x in tree_pruned.tips()]) - set(obs)
        to_delete_rev_set = set(obs) - set([x.name for x in tree_pruned.tips()])
        if len(to_delete_set) > 0 or len(to_delete_rev_set):
            raise ValueError(
                "Pruning the phylogeny failed! There are", len(to_delete_set),
                "features in the phylogeny not present in "
                "the feature table, and",
                len(to_delete_rev_set),
                "features in the feature table not available in the phylogeny!"
                "Both should be 0"
            )
        else:
            print("The phylogeny was pruned successfully!")
    else:
        print("The set of features in the phylogeny and the table "
              "are the same. No feature will be pruned from the tree.")
        tree_pruned = tree
    tree_pruned_out = _1(tree_pruned)

    return tree_pruned_out
