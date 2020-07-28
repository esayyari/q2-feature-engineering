import sys
from qiime2 import Artifact
import pandas as pd
import re
with open(sys.argv[1]) as f:
	all_definition_files = [l.strip() for l in f.readlines()]
	taxonomies_lst = []
	for taxon_file in all_definition_files:
		df_tmp = Artifact.load(taxon_file).view(pd.DataFrame)
		thr = re.sub('clustered_def.', '', taxon_file).replace('.qza', '').replace('./', '')
		df_tmp['threshold.' + thr] = df_tmp.copy('Taxon')
		df_tmp.drop(labels='Taxon', axis=1, inplace=True)
		taxonomies_lst.append(df_tmp)

pd_taxonomies = pd.concat(taxonomies_lst, axis=1)
pd_taxonomies.to_csv(sys.argv[2], sep='\t', index_label='Feature ID')



