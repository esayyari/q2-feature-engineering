# ----------------------------------------------------------------------------
# Copyright (c) 2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Str, Int, Bool, Float, Citations, MetadataColumn, Metadata, Numeric)
from q2_types.feature_table import (FeatureTable, Frequency)
from q2_types.tree import Phylogeny, Rooted
from q2_types.sample_data import SampleData
import re
import ast
import os
from q2_sample_classifier._type import (ClassifierPredictions, Probabilities)

from ._tada import tada

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
                           'original_meta': "The metadata file path corresponding to original samples. "
                                            "This is required when using TADA for balancing. Default is None",
                           'augmented_meta': "The metadata file path corresponding to generated samples. "
                                             "This is required when using TADA for balancing. Default is None"}
_parameters = {'seed_num': Int,
               'meta_data': MetadataColumn[Numeric],
               'xgen': Int,
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
               'original_meta': Str,
               'augmented_meta': Str
               }

_inputs = {'phylogeny': Phylogeny[Rooted],
           'otu_table': FeatureTable[Frequency]}

_input_descriptions = {"phylogeny": "Phylogeny file in newick format",
                       "otu_table": "The count table. This should be in Qiime2 FeatureTable artifact"}

_outputs = [('orig_biom', FeatureTable[Frequency]),
            ('augmented_biom', FeatureTable[Frequency])]

_output_descriptions = {'orig_biom': "Original samples stored in a biom table",
                        'augmented_biom': "Generated samples stored in a biom table"}

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
