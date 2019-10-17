from qiime2.plugin import (Str, Int, Bool, Float, Metadata)
from q2_types.tree import NewickFormat
from qiime2 import CategoricalMetadataColumn
from .TADA_utils import *
import pandas as pd
from qiime2.plugin import TextFileFormat
import biom
import tempfile
from .datasets import make_data_frame
import shutil

def _sort_metada(targets_metadata, biom_table):
    targets = targets_metadata.to_dataframe()

    # filter features and targest so samples match
    index = set(targets.index)
    index = [ix for ix in biom_table.ids('sample') if ix in index]
    targets = targets.loc[index]
    feature_data = biom_table.filter(index, inplace=False)
    return targets, feature_data

def _read_inputs(biom_table: biom.Table, phylogeny_fp: Str, meta_data: CategoricalMetadataColumn = None):
    if meta_data:
        generate_strategy = "balancing"
        meta, biom_table = _sort_metada(meta_data, biom_table)
        y = meta.iloc[:, 0]
        samples = meta.index
    else:
        generate_strategy = "augmentation"
        y = pd.Series(data=np.ones((len(biom_table.ids('sample')), 1)), index=biom_table.ids('sample'))
        samples = biom_table.ids('sample')

    _table = biom_table.sort_order(axis='sample', order=samples)

    _tree = dendropy.Tree.get(path=phylogeny_fp, preserve_underscores=True, schema="newick", rooting='default-rooted')

    if sum(samples != _table.ids('sample')) > 0:
        raise ValueError("The samples IDs in meta data and biom table are not the same! The difference is:",
                         set(samples) - set(_table.ids('sample')), "Please double check.")

    return _table, y, _tree, generate_strategy



def tada(phylogeny: NewickFormat, otu_table: biom.Table, meta_data: CategoricalMetadataColumn = None,
         seed_num: Int = 0, xgen: Int = 0, n_beta: Int = 1, n_binom: Int = 5, var_method: Str = 'br_penalized',
         stat_method: Str = 'binom', prior_weight: Float = 0, coef: Float = 200, exponent: Float = 0.5,
         pseudo_branch_length: Float = 1e-6, pseudo_cnt: Float = 5,
         normalized: Bool = False, output_log_fp: Str = None) -> (biom.Table, biom.Table, pd.Series, pd.Series):
    _table, y, _phylogeny, generate_strategy = _read_inputs(biom_table=otu_table,
                                                            phylogeny_fp=str(phylogeny),
                                                            meta_data=meta_data)
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
        orig_pd, augm_pd = make_data_frame(orig_biom, augm_biom, orig_labels, augm_labels)

    if output_log_fp is not None:
        shutil.copyfile(sG.log_fp, output_log_fp)

    return orig_biom, augm_biom, orig_pd, augm_pd
