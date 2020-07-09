# ----------------------------------------------------------------------------
# Copyright (c) 2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Str, Int, Bool, Float, Citations,
                           MetadataColumn, Numeric, Categorical)
from q2_types.feature_table import (FeatureTable, Frequency)
from q2_types.tree import Phylogeny, Rooted
import re
import ast
import os
from ._tada import tada, prune_features_from_phylogeny
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from ._tada._data_preprocess import cluster_features, reorder_feature_table
from ._smote.ML_over_sampling import synthetic_over_sampling
from ._TreeCluster._treeCluster import tree_cluster
from ._smote.ML_under_sampling import synthetic_under_sampling
from ._CountVectors._countVectors import count_vectors

citations = Citations.load('citations.bib', package='q2_feature_engineering')

_version_re = re.compile(r'__version__\s+=\s+(.*)')

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, '__init__.py'), 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    __version__ = str(ast.literal_eval(hit))

plugin = Plugin(
    name='feature-engineering',
    version=__version__,
    website='https://github.com/tada-alg/TADA',
    package='q2_feature_engineering',
    short_description=('This is a QIIME 2 plugin for phylogenetic augmentation of microbiome samples to '
                       'enhance phenotype classification'),
    citation_text=("Erfan Sayyari, Ban Kawas, Siavash Mirarab, TADA: phylogenetic augmentation of "
                   "microbiome samples enhances phenotype classification, Bioinformatics, "
                   "Volume 35, Issue 14, July 2019, Pages i31–i40, "
                   "https://doi.org/10.1093/bioinformatics/btz394")
)

_parameter_descriptions = {"seed_num": "Seed number. The default value is 0.",
                           "meta_data": "Specifies the generating strategy for either "
                                        "balancing or data augmentation without balancing. "
                                        "If TADA is used for augmentation, this shouldn't "
                                        "be passed. Otherwise, pass a meta data file (in The "
                                        "first column should be samples, and second column "
                                        "should be class labels.",
                           "xgen": "Folds of sample generation for balancing. If TADA is used for only "
                                   "balancing (no extra augmentation afterwards), 0 should be passed. "
                                   "TADA eventually generates new samples until "
                                   "all classes have [xgen+1] * [maximum class size] samples. The "
                                   "default value is 0",
                           "sampling_strategy": "This defines the sampling strategy, the default is None, and will "
                                                "balance the distribution of labels (only for balancing). "
                                                "Alternatively, you can pass a TSV file, where the first column "
                                                "specifies the class label (header column label should be #GroupID), "
                                                "and the second column (with column label 'size'), specifies the size "
                                                "of that class. You should use the same class labels as you have "
                                                "in your metadata column",
                           "n_beta": "The number of draws from the beta distribution. For augmentation, "
                                     "TADA will generate [n_binom]*[n_beta] samples per each sample. "
                                     "The default value is 1.",
                           "n_binom": "The number of draws from binomial distribution. For augmentation,"
                                      " TADA will generate [n_binom]*[n_beta] samples per each sample. "
                                      "The default value is 5",
                           "var_method": "Defines how to introduce the variation. Options are br_penalized"
                                         " and class. The br_penalized can be used with a monotonically "
                                         "increasing function of branch length to define the variation. "
                                         "The class options can be used to use estimate the variation from "
                                         "training data. We suggest using br_penalized (default).",
                           "stat_method": "The generative model. Options are binom or beta_binom, and the "
                                          "default option is binom.",
                           "prior_weight": "The class conditional probability weight. The default is 0.",
                           "coef": "The penalty factor in the calculation of nu. The default value is 200. "
                                   "This affects the amount of variation.",
                           "exponent": "The exponent in the calculation of nu. The default value is "
                                       "0.5. This affects the amount of variation.",
                           "pseudo_branch_length": "A pesudo small branch length will be added to all "
                                                   "branches to avoid zero branch length estimate problem. "
                                                   "The default value is 1e-6.",
                           "pseudo_cnt": "Pseudo count to avoid zero count problem. The default value "
                                         "is adding 5, meaning we add 5/#leaves to each feature value.",
                           "normalized": "If set to 1, the OTU counts will be normalized to add up to "
                                         "one. The default option is 0.",
                           "output_log_fp": "If you want to write the log file of TADA, please specify a path. "
                                            "The default option is None.",
                           'concatenate_meta': "The metadata file path corresponding to the concatenation of "
                                               "generated and original samples. This is required when using TADA for"
                                               "balancing. Default is None",
                           'original_table': "Original samples stored in a biom table. Default is None.",
                           'augmented_table': "Generated samples stored in a biom table. Default is None.",
                           }

_parameters = {'seed_num': Int,
               'meta_data': MetadataColumn[Numeric|Categorical],
               'xgen': Int,
               'sampling_strategy': Str,
               'n_beta': Int,
               'n_binom': Int,
               'var_method': Str,
               'stat_method': Str,
               'prior_weight': Float,
               'coef': Float,
               'exponent': Float,
               'pseudo_branch_length': Float,
               'pseudo_cnt': Int,
               'normalized': Bool,
               'output_log_fp': Str,
               'augmented_table': Str,
               'original_table': Str,
               'concatenate_meta': Str,
               }

_inputs = {'phylogeny': Phylogeny[Rooted],
           'table': FeatureTable[Frequency]}

_input_descriptions = {"phylogeny": "Phylogeny file in newick format",
                       "table": "The count table. This should be in Qiime2 FeatureTable artifact"}

_outputs = [('pruned_phylogeny', Phylogeny[Rooted]),
            ('concatenated_biom', FeatureTable[Frequency])]

_output_descriptions = {'concatenated_biom': "The Concatenation of generated and original biom tables",
                        'pruned_phylogeny': "The pruned input phylogeny based on the training features."}

plugin.methods.register_function(
    function=tada,
    inputs=_inputs,
    parameters=_parameters,
    outputs=_outputs,
    name='Generate microbiome samples using the phylogeny structure.',
    description='Generate microbiome samples with respect to the phylogeny structure.',
    input_descriptions=_input_descriptions,
    parameter_descriptions=_parameter_descriptions,
    output_descriptions=_output_descriptions,
    citations=[citations['TADA']],
    deprecated=False
)

_parameters = {'thr': Float,
               'threads': Int,
               'output_log_file': Str}

_inputs = {'query_table': FeatureTable[Frequency],
           'closed_reference_table': FeatureTable[Frequency],
           'query_sequences': FeatureData[Sequence],
           'reference_sequences': FeatureData[Sequence]}

_input_descriptions = {
    'query_table': 'Dereplicated feature table corresponding to the query sequences',
    'closed_reference_table': 'Feature table corresponding to the previously',
    'query_sequences': 'Query dereplicated fragments',
    'reference_sequences': 'The sequences to use as cluster centroids'
}

_outputs = [('clustered_table', FeatureTable[Frequency]),
            ('clustered_sequences', FeatureData[Sequence]),
            ('unmatched_sequences', FeatureData[Sequence])]

_output_descriptions = {'clustered_table': "The table following clustering of features.",
                        'clustered_sequences': "The sequences representing clustered features, "
                                               "relabeled by the reference IDs.",
                        'unmatched_sequences': "The sequences which failed to match any reference "
                                               "sequences. This output maps to vsearch's"}

_parameter_descriptions = {'thr': 'The similarity threshold for clustering.',
                           'threads': 'The number of threads to be used for clustering. '
                                      'The default is 0 (total number of available threads)',
                           'output_log_file': 'The path to save the log file of the code. Default is None'}

plugin.methods.register_function(
    function=cluster_features,
    inputs=_inputs,
    parameters=_parameters,
    outputs=_outputs,
    name='Does closed reference OTU clustering using a limited reference.',
    description="This tool first prunes the reference database using a clustered feature table, "
                "and then run closed reference OTU picking for the query sequences using the pruned reference database.",
    input_descriptions=_input_descriptions,
    parameter_descriptions=_parameter_descriptions,
    output_descriptions=_output_descriptions,
    citations=[citations['rognes2016vsearch']],
    deprecated=False
)

plugin.methods.register_function(
    function=reorder_feature_table,
    inputs={'query_table': FeatureTable[Frequency],
            'reference_table': FeatureTable[Frequency]},
    outputs=[('reordered_table', FeatureTable[Frequency])],
    parameters={},
    parameter_descriptions={},
    input_descriptions={'query_table': 'Input frequency table to be reordered.',
                        'reference_table': 'Reference table that specifies the order of features'},
    output_descriptions={'reordered_table': "The reordered table."},
    name='This is a plugin for reordering a feature table based on another table.',
    description='This is a plugin for reordering a feature table based on another table.',
    citations=[],
    deprecated=False
)


_inputs = {'table': FeatureTable[Frequency]}
_parameters = {'metadata': MetadataColumn[Numeric|Categorical],
               'method': Str,
               'k_neighbors': Int,
               'n_jobs': Int,
               'sampling_strategy': Str,
               'random_state': Int,
               'output_log_fp': Str,
               'concatenate_meta_fp': Str}

_outputs = [('over_sampled_table', FeatureTable[Frequency])]

_output_descriptions = {'over_sampled_table': "The Concatenation of generated and original biom tables"}

_parameter_descriptions = {'metadata': "Specifies the generating strategy for either "
                                       "balancing or data augmentation without balancing. "
                                       "If TADA is used for augmentation, this shouldn't "
                                       "be passed. Otherwise, pass a meta data file (in The "
                                       "first column should be samples, and second column "
                                       "should be class labels.",
                           'method': "The over sampling method to be used. The default is SMOTE and options are:"
                                     "ADASYN, RandomOverSampler, and SMOTE."
                                     "Using 'RandomOverSampler', and specifying the correct sampling strategy, "
                                     "we can also do down-sampling.",
                           'k_neighbors': "number of nearest neighbours to used to construct synthetic samples. "
                                          "The default value is 5",
                           'n_jobs': 'The number of threads to be used for clustering. Default is 1.',
                           'sampling_strategy': "This defines the sampling strategy, the default is 'auto', "
                                                "and will balance the distribution of labels. "
                                                "Alternatively, you can pass a TSV file, where the first "
                                                "column specifies the class label "
                                                "(header column label should be #SampleID), and the second column "
                                                "(with column label size), specifies the size of that class.",
                           'random_state': 'Seed number. The default value is 42.',
                           "output_log_fp": "If you want to write the log file of TADA, please specify a path. "
                                            "The default option is None.",
                           'concatenate_meta_fp': "The metadata file path corresponding to the concatenation of"
                           }
_input_descriptions = {'table': 'The feature table that will be used for over sampling.'}

plugin.methods.register_function(
    function=synthetic_over_sampling,
    inputs=_inputs,
    outputs=_outputs,
    parameters=_parameters,
    parameter_descriptions=_parameter_descriptions,
    input_descriptions=_input_descriptions,
    output_descriptions=_output_descriptions,
    name="The plugin for balancing class labels by over sampling the minority class.",
    description="The plugin for balancing class labels by over sampling the minority class.",
    citations=[citations['imblearn'], citations['ADASYN'], citations['chawla2002smote']],
    deprecated=False
)

_inputs = {'table': FeatureTable[Frequency]}
_parameters = {'metadata': MetadataColumn[Numeric|Categorical],
               'method': Str,
               'voting': Str,
               'n_jobs': Int,
               'sampling_strategy': Str,
               'random_state': Int,
               'output_log_fp': Str,
               'concatenate_meta_fp': Str}

_outputs = [('under_sampled_table', FeatureTable[Frequency])]

_output_descriptions = {'under_sampled_table': "The Concatenation of generated and original biom tables"}

_parameter_descriptions = {'metadata': "Specifies the generating strategy for either "
                                       "balancing or data augmentation without balancing. "
                                       "If TADA is used for augmentation, this shouldn't "
                                       "be passed. Otherwise, pass a meta data file (in The "
                                       "first column should be samples, and second column "
                                       "should be class labels.",
                           'method': "The over sampling method to be used. The default is RandomUnderSampler and "
                                     "option is RandomUnderSampler.",
                           'voting': " Voting strategy to generate the new samples. Options are soft, hard and "
                                     "the default is auto.",
                           'n_jobs': 'The number of threads to be used for clustering. Default is 1.',
                           'sampling_strategy': "This defines the sampling strategy, the default is 'auto', "
                                                "and will balance the distribution of labels. "
                                                "Alternatively, you can pass a TSV file, where the first "
                                                "column specifies the class label "
                                                "(header column label should be #SampleID), and the second column "
                                                "(with column label size), specifies the size of that class.",
                           'random_state': 'Seed number. The default value is 42.',
                           "output_log_fp": "If you want to write the log file of TADA, please specify a path. "
                                            "The default option is None.",
                           'concatenate_meta_fp': "The metadata file path corresponding to the concatenation of"
                           }
_input_descriptions = {'table': 'The feature table that will be used for under sampling.'}

plugin.methods.register_function(
    function=synthetic_under_sampling,
    inputs=_inputs,
    outputs=_outputs,
    parameters=_parameters,
    parameter_descriptions=_parameter_descriptions,
    input_descriptions=_input_descriptions,
    output_descriptions=_output_descriptions,
    name="The plugin for balancing class labels by under sampling the over represented class.",
    description="The plugin for balancing class labels by under sampling the over represented class.",
    citations=[citations['imblearn']],
    deprecated=False
)

_inputs = {'table': FeatureTable[Frequency],
           'phylogeny': Phylogeny[Rooted]}

_outputs = [('pruned_phylogeny', Phylogeny[Rooted])]
_parameters = {}
_parameter_descriptions = {}
_output_descriptions = {'pruned_phylogeny': 'Pruned phylogeny'}
_input_descriptions = {'table': 'Feature table',
                       'phylogeny': "The phylogeny corresponding to the set of features available in the table. "
                                    "There can be features not present in the table that will be removed from the tree"}
plugin.methods.register_function(
    function=prune_features_from_phylogeny,
    inputs=_inputs,
    outputs=_outputs,
    parameters=_parameters,
    parameter_descriptions=_parameter_descriptions,
    input_descriptions=_input_descriptions,
    output_descriptions=_output_descriptions,
    name="Tool to prune a phylogeny given the feature table",
    description="Tool to prune a phylogeny using the features present in the feature table. ",
    citations=[],
    deprecated=False
)

_inputs = {'phylogeny_fp': Phylogeny[Rooted],
           'table': FeatureTable[Frequency]}

_input_descriptions = {'phylogeny_fp': 'The tree with inserted feature data',
                       'table': 'The count table'}


_outputs = [('clustered_phylogeny', Phylogeny[Rooted]),
            ('clusters_definitions', FeatureData[Taxonomy]),
            ('clustered_table', FeatureTable[Frequency])]


_output_descriptions = {'clustered_phylogeny': 'The tree with leaves defined based on the cluster definitions',
                        'clusters_definitions': 'Cluster definition',
                        'clustered_table': 'The count table defined based on the cluster definitions'}

_parameters = {'method': Str,
               'threshold': Float,
               'threshold_free': Str}

_parameter_descriptions = {'method': 'The method used for clustering leaves. Default is max_clade',
                           'threshold': 'The threshold used for clustering leaves. Default is 0',
                           'threshold_free': 'To use TreeCluster in the threshold free mode set it as argmax_clusters.'
                                             'The default is None'}


plugin.methods.register_function(
    function=tree_cluster,
    inputs=_inputs,
    outputs=_outputs,
    parameters=_parameters,
    parameter_descriptions=_parameter_descriptions,
    input_descriptions=_input_descriptions,
    output_descriptions=_output_descriptions,
    name="TreeCluster",
    description="TreeCluster is a tool that, given a tree T (Newick format) and a distance threshold t, "
                "finds the minimum number of clusters of the leaves of T such that some user-specified "
                "constraint is met in each cluster. ",
    citations=[citations['TreeCluster']],
    deprecated=False
)


_inputs = {'phylogeny': Phylogeny[Rooted],
           'table': FeatureTable[Frequency]}

_input_descriptions = {'phylogeny': 'The tree with inserted feature data',
                       'table': 'The count table'}


_outputs = [('count_vector_table', FeatureTable[Frequency])]


_output_descriptions = {'count_vector_table': 'The count table defined based '
                                              'on the cluster definitions'}

_parameters = {'method': Str}

_parameter_descriptions = {'method': 'The method used for calculating '
                                     'count vectors. Options are '
                                     'weighted_unifrac or unweighted_unifrac. '
                                     'The default is weighted_unifrac'}


plugin.methods.register_function(
    function=count_vectors,
    inputs=_inputs,
    outputs=_outputs,
    parameters=_parameters,
    parameter_descriptions=_parameter_descriptions,
    input_descriptions=_input_descriptions,
    output_descriptions=_output_descriptions,
    name="CountVectors",
    description="This plugin is to extract weighted and unweighted count "
                "vectors using a phylogenetic tree ",
    citations=[citations['FastUniFrac']],
    deprecated=False
)


