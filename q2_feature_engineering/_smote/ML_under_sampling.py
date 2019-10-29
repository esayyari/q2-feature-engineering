from imblearn import under_sampling
from qiime2.plugin import (Str, Int)
import biom
from q2_feature_engineering._tada.logger import LOG
from qiime2 import NumericMetadataColumn
import numpy as np
import pandas as pd
import qiime2
import tempfile
import shutil

dispatcher = {'RandomUnderSampler': under_sampling.RandomUnderSampler}


def _sort_metada(targets_metadata, biom_table):
    targets = targets_metadata.to_dataframe()

    # filter features and targest so samples match
    index = set(targets.index)
    index = [ix for ix in biom_table.ids('sample') if ix in index]
    targets = targets.loc[index]
    feature_data = biom_table.filter(index, inplace=False)
    return targets, feature_data


def _read_inputs(biom_table: biom.Table, meta_data: NumericMetadataColumn = None):
    if meta_data:
        meta, biom_table = _sort_metada(meta_data, biom_table)
        y = meta.iloc[:, 0]
        samples = meta.index
    else:
        samples = biom_table.ids('sample')
        y = pd.DataFrame(data=np.asarray(np.ones((len(samples), 1))).ravel(), index=samples)

    _table = biom_table.sort_order(axis='sample', order=samples)

    if np.sum(samples != _table.ids('sample')) > 0:
        raise ValueError("The samples IDs in meta data and biom table are not the same! The difference is:",
                         set(samples) - set(_table.ids('sample')), "Please double check.")

    return _table, y


def synthetic_under_sampling(table: biom.Table, metadata: NumericMetadataColumn,
                             concatenate_meta_fp: Str, method: Str = 'RandomUnderSampler',
                             voting: Str = 'auto', n_jobs: Int = 1,
                             sampling_strategy: Str = 'auto',
                             random_state: Int = 42, output_log_fp: Str = None) -> biom.Table:
    log_fp = tempfile.mktemp()
    print("The log file will be writen into", log_fp)
    if log_fp:
        logger_ins = LOG(log_fp=log_fp).get_logger('synthetic_over_sampling')
        logger_ins.info("The parameters used for oversampling are")
        logger_ins.info('voting (will be used with ClusterCentroids only):', voting)
        logger_ins.info('Sampling method:', method)
        logger_ins.info('Output log file path:', log_fp)
        logger_ins.info('sampling_strategy:', sampling_strategy)
        logger_ins.info('n_jobs:', n_jobs)
        logger_ins.info('random_state:', random_state)

    cls = dispatcher[method]
    if method != 'RandomUnderSampler':
        table.norm(inplace=True)
        if log_fp:
            logger_ins.info("The input table is normalized before using it for oversampling")
    sorted_table, sorted_metadata = _read_inputs(table, meta_data=metadata)
    matrix_data = sorted_table.matrix_data.transpose().todense()
    if method not in dispatcher:
        raise ValueError(
            'The optional methods for over sampling are', dispatcher.keys(), "instead it received", method
        )
    if method == 'RandomUnderSampler':
        under_sampling_cls = cls(sampling_strategy=sampling_strategy, random_state=random_state, replacement=False)
    else:
        raise NotImplementedError("Method", method, "is not implemented yet")
    X_resampled, y_resampled = under_sampling_cls.fit_resample(matrix_data, sorted_metadata)
    if log_fp:
        logger_ins.info("The under-sampling finished successfully!")
        logger_ins.info("Overall, the size of data is", len(X_resampled))
    if method == 'RandomUnderSampler':
        dummy_samples = under_sampling_cls.sample_indices_
    else:
        raise NotImplementedError("Method", method, "is not implemented yet")
    under_sampling_dummy = sorted_table.filter(ids_to_keep=dummy_samples, inplace=False)
    under_sampling_dummy = under_sampling_dummy.sort_order(order=dummy_samples, axis='sample')
    if method == "RandomUnderSampler" and np.sum(under_sampling_dummy.matrix_data-X_resampled) != 0:
        raise ValueError("The undersampling changed the matrix data")

    undersampled_table = biom.Table(X_resampled.transpose(), observation_ids=sorted_table.ids('observation'),
                                    sample_ids=dummy_samples)
    undersampled_metadata = pd.DataFrame(index=dummy_samples, data=y_resampled)
    undersampled_metadata.index.names = ['#SampleID']
    undersampled_metadata.columns = ['label']
    undersampled_meta = qiime2.Metadata(undersampled_metadata)
    undersampled_meta.save(concatenate_meta_fp)

    if log_fp:
        shutil.copy(log_fp, output_log_fp)

    return undersampled_table
