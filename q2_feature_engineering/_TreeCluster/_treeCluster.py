import numpy as np
import biom
import pandas as pd


def initialize_mapping_matrix(mapping, table):
    m = max(mapping.ClusterNumber)
    n = len(table.ids('observation'))
    k = max(mapping.ClusterNumber)
    map_matrix = np.zeros((m, n))
    return map_matrix


def make_mapping_matrix(mapping, map_matrix):
    i = 0
    j = 0
    obs = [""] * len(map_matrix)
    singletones = []
    for idx, row in mapping.iterrows():
        if row['ClusterNumber'] == -1:
            singletones.append(i)
            j += 1
            obs.append(idx)
        else:
            map_matrix[row['ClusterNumber']-1][i] = 1
            obs[row['ClusterNumber']-1] = 'cluster' + str(row['ClusterNumber'])
        i += 1
    singletones = np.asarray(singletones)
    return map_matrix, singletones, obs


def check_mapping(final_data, table, new_data, mapping, map_matrix):
    flag = False
    check = []

    if not final_data.shape[0] == new_data.shape[0] + (mapping.ClusterNumber == -1).sum():
        print("The number of features in the mapped and original tables are inconsistent!")
        print(final_data.shape[0], new_data.shape[0] + (mapping.ClusterNumber == -1).sum())
        flag = True
    if not final_data.sum() == table.matrix_data.sum():
        print("The sum of all faetures in the original table and mapped one are not the same!")
        print(final_data.sum(), table.matrix_data.sum())
        flag = True
    for i in range(new_data.shape[0]):
        check.append(set(mapping.iloc[np.where(map_matrix[i,:])].ClusterNumber.values) == {i + 1})
    if np.sum(check) != len(check):
        print("There is an inconsistency in the mapped features")
        flag = True
    if not flag:
        print("The feature mapping finished successfully! ")
    return flag


def map_features(mapping: pd.DataFrame, table: biom.Table) -> biom.Table:
    mapping = mapping.reset_index().set_index('SequenceName')
    mapping = mapping.loc[table.ids('observation')].reindex(table.ids('observation'))
    map_matrix = initialize_mapping_matrix(mapping, table)
    map_matrix, singletones, obs = make_mapping_matrix(mapping, map_matrix)
    new_data = map_matrix * table.matrix_data
    final_data = np.row_stack([new_data, np.asarray(table.matrix_data[singletones, :].todense())])
    mapped_biom_table = biom.Table(data=final_data, observation_ids=obs, sample_ids=table.ids())
    if check_mapping(final_data, table, new_data, mapping, map_matrix):
        raise ValueError("something went wrong!")
    return mapped_biom_table

