# ----------------------------------------------------------------------------
# Copyright (c) 2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.tree import NewickFormat
from qiime2.plugin import (Int, Bool)
from qiime2.plugins import feature_table
import biom
import pandas as pd
import dendropy


def _update_seq_feature_ids(map_dict: dict, seqs_df: pd.Series) -> pd.Series:
    seqs_dict = seqs_df.to_dict()
    mapped_seqs_dict = dict()
    for sp, seq in seqs_dict.items():
        mapped_sp = map_dict[sp]
        mapped_seqs_dict[mapped_sp] = seq
    return pd.Series(mapped_seqs_dict)


def _create_inv_seq_mapping(seqs_df: pd.Series):
    mapping = dict()
    for sp, seq in seqs_df.iteritems():
        mapping[seq] = sp
    return mapping


def _update_table_feature_ids(mapping: dict, table: biom.Table) -> biom.Table:
    return table.update_ids(mapping, axis='observation', inplace=False)


def _update_table_sample_ids(mapping: dict, table: biom.Table) -> biom.Table:
    return table.update_ids(mapping, axis='sample', inplace=False)


def _prune_phylo(phylogeny: dendropy.Tree, species_set: set) -> dendropy.Tree:
    to_delete_set = set()
    for sp in species_set:
        to_delete_set.add(phylogeny.find(sp))
    phylogeny.remove_deleted(lambda x: x in to_delete_set)
    phylogeny.prune()
    return phylogeny


def _prune_table(table: biom.Table, min_samples: Int = 1) -> biom.Table:
    filtered_table = feature_table.methods.filter_features(table, min_samples=min_samples)
    return filtered_table.filtered_table


def _prune_seqs(seqs: pd.Series, table: biom.Table, verbose: Bool = False) -> pd.Series:
    if verbose:
        print("The total number of sequences is", seqs.size)
    filtered_seqs = feature_table.methods.filter_seqs(seqs, table=table)
    filtered_seqs = filtered_seqs.filtered_data
    if verbose:
        print("The number of sequences after applying the filter is now", filtered_seqs.view(pd.Series).size)
    return filtered_seqs


def _select_samples(table: biom.Table, metadata: pd.DataFrame, verbose: Bool = False) -> biom.Table:
    if verbose:
        print("The total number of samples is", len(table.view(biom.Table).ids('sample')))
        print("The total number of features is", len(table.view(biom.Table).ids('observation')))
    filtered_table = feature_table.methods.filter_samples(table=table, metadata=metadata)
    filtered_table = filtered_table.filtered_table
    if verbose:
        print("The number of samples is now", len(filtered_table.view(biom.Table).ids('sample')))
    filtered_table = feature_table.methods.filter_features(table=filtered_table, min_samples=1)
    filtered_table = filtered_table.filtered_table
    if verbose:
        print("The number of features after removing samples is ",
              len(filtered_table.view(biom.Table).ids('observation')))
    return filtered_table



