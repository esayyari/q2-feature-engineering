# ----------------------------------------------------------------------------
# Copyright (c) 2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


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
from .logger import LOG


def _sort_metada(targets_metadata, biom_table):
    targets = targets_metadata.to_dataframe()

    # filter features and targest so samples match
    index = set(targets.index)
    index = [ix for ix in biom_table.ids('sample') if ix in index]
    targets = targets.loc[index]
    feature_data = biom_table.filter(index, inplace=False)
    return targets, feature_data


def _read_inputs(biom_table: biom.Table, phylogeny_fp: Str, meta_data: NumericMetadataColumn = None):
    if meta_data:
        generate_strategy = "balancing"
        meta, biom_table = _sort_metada(meta_data, biom_table)
        y = meta.iloc[:, 0]
        samples = meta.index
    else:
        generate_strategy = "augmentation"
        y = pd.Series(data=np.ones((len(biom_table.ids('sample')),)), index=biom_table.ids('sample'))
        samples = biom_table.ids('sample')

    _table = biom_table.sort_order(axis='sample', order=samples)

    _tree = dendropy.Tree.get(path=phylogeny_fp, preserve_underscores=True, schema="newick", rooting='default-rooted')

    if sum(samples != _table.ids('sample')) > 0:
        raise ValueError("The samples IDs in meta data and biom table are not the same! The difference is:",
                         set(samples) - set(_table.ids('sample')), "Please double check.")

    return _table, y, _tree, generate_strategy


def tada(phylogeny: NewickFormat, otu_table: biom.Table, meta_data: NumericMetadataColumn = None,
         seed_num: Int = 0, xgen: Int = 0, n_beta: Int = 1, n_binom: Int = 5, var_method: Str = 'br_penalized',
         stat_method: Str = 'binom', prior_weight: Float = 0, coef: Float = 200, exponent: Float = 0.5,
         pseudo_branch_length: Float = 1e-6, pseudo_cnt: Float = 5,
         normalized: Bool = False, output_log_fp: Str = None,
         augmented_meta: Metadata = None, original_meta: Metadata = None,
         concatenate_meta: Metadata = None) -> (biom.Table, biom.Table, biom.Table):
    _table, y, _phylogeny, generate_strategy = _read_inputs(biom_table=otu_table,
                                                            phylogeny_fp=str(phylogeny),
                                                            meta_data=meta_data)
    if generate_strategy is 'balancing' and (augmented_meta is None or original_meta is None):
        raise ValueError(
            "Expected a path to write out the generated and original labels and metadata!"
        )
    with tempfile.TemporaryDirectory() as tmp:
        sG = SampleGenerator(seed_num=seed_num, logger=None, generate_strategy=generate_strategy, tmp_dir=tmp,
                             xgen=xgen, n_beta=n_beta, n_binom=n_binom, var_method=var_method, stat_method=stat_method,
                             prior_weight=prior_weight,
                             coef=coef, exponent=exponent, pseudo=pseudo_branch_length,
                             pseudo_cnt=pseudo_cnt, normalized=normalized)
        orig_biom, orig_labels, augm_biom, augm_labels = sG.fit_transform(table=_table, y=y, tree=_phylogeny)
        if np.sum(orig_biom.matrix_data - otu_table.matrix_data) > 1e-20:
            raise ValueError(
                "The original biom table doesn't match the output of generator function! Please double check")
        if generate_strategy is 'balancing':
            orig_pd, augm_pd = make_data_frame(orig_biom, augm_biom, orig_labels, augm_labels)
            if not os.path.exists(os.path.dirname(str(original_meta))):
                os.mkdir(os.path.dirname(str(original_meta)))
            orig_meta = qiime2.Metadata(orig_pd)
            augm_meta = qiime2.Metadata(augm_pd)
            concat_pd = pd.concat([orig_pd, augm_pd])
            concat_meta = qiime2.Metadata(concat_pd)
            orig_meta.save(original_meta)
            augm_meta.save(augmented_meta)
            concat_meta.save(concatenate_meta)

        if output_log_fp is not None:
            if not os.path.exists(os.path.dirname(output_log_fp)):
                os.mkdir(os.path.dirname(output_log_fp))
            shutil.copyfile(sG.log_fp, output_log_fp)
        if np.sum(orig_biom.ids('observation') == augm_biom.ids('observation')) != len(orig_biom.ids('observation')):
            raise ValueError(
                "The order of features in original and augmented data is different."
                "Please make sure that your phylogeny doesn't have extra features"
            )
        concat_biom = orig_biom
        concat_biom = concat_biom.merge(augm_biom)

    return orig_biom, augm_biom, concat_biom


def prune_features_from_phylogeny(table: biom.Table, phylogeny: NewickFormat, out_log_fp: Str = None) -> NewickFormat:
    log_fp = tempfile.mktemp()
    logger_ins = LOG(log_fp).get_logger('prune_phylogeny')
    tree = dendropy.Tree.get(path=str(phylogeny), preserve_underscores=True, schema="newick", rooting='default-rooted')
    obs = table.ids('observation')
    to_delete_set = set([x.taxon.label for x in tree.leaf_nodes()]) - set(obs)
    if len(set(obs) - set([x.taxon.label for x in tree.leaf_nodes()])) > 0:
        raise ValueError(
            "There are",  len(set(obs) - set([x.taxon.label for x in tree.leaf_nodes()])),
            "features in the feature table not present in the phylogeny! Please check your tree"
        )
    if len(to_delete_set) > 0:
        logger_ins.info("The set of features in the phylogeny and the table are not the same.",
                        len(to_delete_set), "features will be pruned from the tree.")
        tree.prune_taxa_with_labels(to_delete_set)
    else:
        logger_ins.info("The set of features in the phylogeny and the table are the same. "
                        "No feature will be pruned from the tree.")
    if log_fp:
        shutil.copy(log_fp, out_log_fp)
    tree_out = NewickFormat()
    tmp_out = tempfile.mktemp()
    tree.write(path=tmp_out, schema="newick", suppress_rooting=True)
    shutil.copy(tmp_out, str(tree_out))
    return tree_out
