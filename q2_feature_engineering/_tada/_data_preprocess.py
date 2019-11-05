import biom
import pandas as pd
import numpy as np
from q2_vsearch._cluster_features import cluster_features_closed_reference
from q2_types.feature_data import DNAFASTAFormat
from qiime2 import Artifact
from q2_types.feature_data._transformer import _16, _15
import tempfile
from q2_feature_engineering._tada.logger import LOG
import shutil


def get_reference_seqs_from_ids(table: biom.Table, reference_seqs_pd: pd.Series) -> DNAFASTAFormat:
    output_references = pd.Series()
    for obs in table.ids('observation'):
        seq = reference_seqs_pd[obs]
        output_references[obs] = seq
    output_references_fasta = _16(output_references)
    return output_references_fasta


def cluster_features(query_table: biom.Table, closed_reference_table: biom.Table,
                                         query_sequences: DNAFASTAFormat,
                                         reference_sequences: pd.Series, thr: float = 0.97,
                                         threads: int = 1, output_log_file: str = None) -> (
        biom.Table, DNAFASTAFormat, DNAFASTAFormat):
    reference_sequences_fasta = get_reference_seqs_from_ids(closed_reference_table, reference_sequences)
    results = cluster_features_closed_reference(sequences=query_sequences, table=query_table,
                                                reference_sequences=reference_sequences_fasta,
                                                perc_identity=thr, threads=threads)

    clustered_table_biom = results[0]

    clustered_sequences_pd = Artifact.load(str(results[1])).view(pd.Series)
    unmatched_sequences_pd = Artifact.load(str(results[2])).view(pd.Series)

    with tempfile.mktemp() as tmp_fp:
        logger_ins = LOG(tmp_fp).get_logger('clustering_features')
        logger_ins.info("The number of OTUs in the reference database is", _15(reference_sequences_fasta).size)
        logger_ins.info("The number of unmatched sequence to the reference alignment is", unmatched_sequences_pd.size)
        logger_ins.info("The number of matched sequences to the reference alignment is", clustered_sequences_pd.size)
        logger_ins.info("Before applying clustering, the total number of counts "
                        "in the original feature table was", np.sum(query_table.sum()))
        logger_ins.info("Before applying clustering, the number of non-zero elements"
                        " of the underlying feature table is", query_table.nnz)
        logger_ins.info("After applying clustering, the total number of counts "
                        "in the original feature table was", np.sum(clustered_table_biom.sum()))
        logger_ins.info("After applying clustering, the number of non-zero elements"
                        " of the underlying feature table is", clustered_table_biom.nnz)
        logger_ins.info("The percent of total counts retained is",
                        np.sum(query_table.sum()) / np.sum(clustered_table_biom.sum()) * 100, "%s")

        query_samples = clustered_table_biom.ids('sample')
        closed_reference_features = closed_reference_table.ids('observation')
        clustered_table_biom = closed_reference_table.merge(clustered_table_biom)
        clustered_table_biom.filter(ids_to_keep=query_samples, axis='sample', inplace=True)
        if len(set(closed_reference_features) - set(clustered_table_biom.ids('sample'))) != 0:
            raise ValueError(
                "Merging two tables failed! There are less features in the final table than expected!"
            )
        if output_log_file:
            shutil.copy(tmp_fp, output_log_file)
    return clustered_table_biom, results[1], results[2]


def reorder_feature_table(query_table: biom.Table, reference_table: biom.Table) -> biom.Table:
    query_samples = query_table.ids()
    ref_samples = set(reference_table.ids())
    for sample in query_samples:
        if sample in ref_samples:
            raise ValueError(
                "The sample", sample, "from your reference data found in your query one, while "
                                      "the two tables should be disjoint."
            )
    merged_table = reference_table.merge(query_table)
    merged_table.filter(ids_to_keep=reference_table.ids('observation'))
    merged_table.filter(ids_to_keep=query_samples)
    for sample in merged_table.ids():
        if sample in ref_samples:
            raise ValueError(
                "The sample", sample, "from your reference data found in your query one, while "
                                      "the two tables should be disjoint."
            )

    return merged_table.sort_order(reference_table.ids('observation'))

