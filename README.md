# Extracting and preparing training data using q2-feature-engineerng
 
In this tutorial, we will introduce `q2-feature-engineering` for extracting and manipulating data used for training supervised or unsupervised machine learning methods. For now, we introduce a plugin for data augmentation, TADA (Sayyari et. al., Bioinformatics, 2019). In this tutorial, we will use a subset of American Gut Project (AGP, Daniel McDonald et. al. mSystems, 2018) to generate new data using TADA for classifying Body Mass Index (BMI) categories. 


## Introduction
TADA is a new data augmentation technique for classifying phenotypes based on the microbiome. Our algorithm, TADA, uses training biological data and a statistical generative model to create new samples adding to the existing ones, addressing issues of low-sample-size. In generating new samples, TADA takes into account phylogenetic relationships between microbial species. Adding these synthetic samples to the training set improves the accuracy of downstream classification, especially when the training data have an unbalanced representation of classes.


## Installation
q2-feature-engineering is a Qiime2 plugin, but is not officially integrated with other plugins. This python3 plugin depends on several other python packages including Qiime2 itself. These dependencies are i) [numpy](https://www.numpy.org/) ii) [Dendropy](https://dendropy.org/) iii) [scikit-learn](https://scikit-learn.org/stable/install.html) iv) [biom-format](http://biom-format.org/documentation/biom_format.html) v) [pandas](https://pandas.pydata.org/) vi) [Qiime2](https://qiime2.org/). 

You need to first install Qiime2, following this [instruction](https://docs.qiime2.org/2019.7/install/). We created this tutorial using Core 2019.7 distribution (conda instruction), and so we assume that the conda environment is `qiime2-2019.7`. By installing Qiime2 [numpy](https://www.numpy.org/), [biom-format](http://biom-format.org/documentation/biom_format.html), [pandas](https://pandas.pydata.org/) requirements will be automatically fulfilled. You only need to further install [Dendropy](https://dendropy.org/) and [scikit-learn](https://scikit-learn.org/stable/install.html) 


```
conda install scikit-learn
conda install -c bioconda dendropy
```

Now activate your conda environment using either of your commands (which ever works for your current installation of conda)

```
conda activate qiime2-2019.7
or 
source activate  qiime2-2019.7
```

Next clone this [github](git@github.com:tada-alg/TADA.git) repository somewhere on your machine:

```
cd [Any directory you wish to install TADA]
git clone git@github.com:tada-alg/TADA.git
cd TADA
```

You are now ready to install the package:

```
python setup.py install
```

This will install Qiime2 plugin as well as standalone TADA package for you. You can test your installation using the following commands:

```
# Update Qiime package
qiime dev refresh-cache

# See intputs to TADA
qiime feature-engineering tada --help



Usage: qiime feature-engineering tada [OPTIONS]

  Generate microbiome samples with respect to the phylogeny structure.

Inputs:
  --i-phylogeny ARTIFACT  Phylogeny file in newick format
    Phylogeny[Rooted]                                               [required]
  --i-otu-table ARTIFACT FeatureTable[Frequency]
                          The count table. This should be in Qiime2
                          FeatureTable artifact                     [required]
Parameters:
  --m-meta-data-file METADATA
  --m-meta-data-column COLUMN  MetadataColumn[Numeric]
                          Specifies the generating strategy for either
                          balancing or data augmentation without balancing. If
                          TADA is used for augmentation, this shouldn't be
                          passed. Otherwise, pass a meta data file (in The
                          first column should be samples, and second column
                          should be class labels.                   [optional]
  --p-seed-num INTEGER    Seed number. The default value is 0.    [default: 0]
  --p-xgen INTEGER        Folds of sample generation for balancing. If TADA
                          is used for only balancing (no extra augmentation
                          afterwards), 0 should be passed. TADA eventually
                          generates new samples until all classes have
                          [xgen+1] * [maximum class size] samples. The default
                          value is 0                              [default: 0]
  --p-n-beta INTEGER      The number of draws from the beta distribution. For
                          augmentation, TADA will generate [n_binom]*[n_beta]
                          samples per each sample. The default value is 1.
                                                                  [default: 1]
  --p-n-binom INTEGER     The number of draws from binomial distribution. For
                          augmentation, TADA will generate [n_binom]*[n_beta]
                          samples per each sample. The default value is 5
                                                                  [default: 5]
  --p-var-method TEXT     Defines how to introduce the variation. Options are
                          br_penalized and class. The br_penalized can be used
                          with a monotonically increasing function of branch
                          length to define the variation. The class options
                          can be used to use estimate the variation from
                          training data. We suggest using br_penalized
                          (default).                 [default: 'br_penalized']
  --p-stat-method TEXT    The generative model. Options are binom or
                          beta_binom, and the default option is binom.
                                                            [default: 'binom']
  --p-prior-weight NUMBER The class conditional probability weight. The
                          default is 0.                           [default: 0]
  --p-coef NUMBER         The penalty factor in the calculation of nu. The
                          default value is 200. This affects the amount of
                          variation.                            [default: 200]
  --p-exponent NUMBER     The exponent in the calculation of nu. The default
                          value is 0.5. This affects the amount of variation.
                                                                [default: 0.5]
  --p-pseudo-branch-length NUMBER
                          A pesudo small branch length will be added to all
                          branches to avoid zero branch length estimate
                          problem. The default value is 1e-6. [default: 1e-06]
  --p-pseudo-cnt INTEGER  Pseudo count to avoid zero count problem. The
                          default value is adding 5, meaning we add 5/#leaves
                          to each feature value.                  [default: 5]
  --p-normalized / --p-no-normalized
                          If set to 1, the OTU counts will be normalized to
                          add up to one. The default option is 0.
                                                              [default: False]
  --p-output-log-fp TEXT  If you want to write the log file of TADA, please
                          specify a path. The default option is None.
                                                                    [optional]
  --p-augmented-meta TEXT The metadata file path corresponding to generated
                          samples. This is required when using TADA for
                          balancing. Default is None                [optional]
  --p-original-meta TEXT  The metadata file path corresponding to original
                          samples. This is required when using TADA for
                          balancing. Default is None                [optional]
Outputs:
  --o-orig-biom ARTIFACT FeatureTable[Frequency]
                          Original samples stored in a biom table   [required]
  --o-augmented-biom ARTIFACT FeatureTable[Frequency]
                          Generated samples stored in a biom table  [required]
Miscellaneous:
  --output-dir PATH       Output unspecified results to a directory
  --verbose / --quiet     Display verbose output to stdout and/or stderr
                          during execution of this action. Or silence output
                          if execution is successful (silence is golden).
  --citations             Show citations and exit.
  --help                  Show this message and exit.
```

For the purpose of this tutorial, now please go to `data/reference` folder

```
cd data/reference
tar xzvf test.tar.gz
```

This will create the following folder/files that

```
test/
test/metadata.csv
test/phylogeny.tre
test/feature-table.biom
```

The main inputs to TADA are:

1. A rooted phylogeny __pruned__ to the fragments available in your training dataset (Qiime2 Phylogeny[Rooted] artifact, please see below for more information). For microbiome data, you could use this [pipline](https://github.com/qiime2/q2-fragment-insertion) to get the tree.
2. A biom table contains your sOTU/OTU __counts__ (it doesn't work with normalized counts, i.e. frequencies) (Qiime2 FeatureTable[Frequency] artifact, please see below for more information).
3. If Beta-Binomial generative model will be used for augmentation, you need to pass a metadata file (comma seperated) with a header line. The first column corresponds to subject IDs (please use _#SampleID_ as the header), followed by numerical (integer) columns (use any desired labels with no tabs for column labels). These columns can be anything from categories of phenotypes (converted to integers) to unsupervised labels. (imported as SampleData[ClassifierPredictions], please see below for more information)
4. If you are passing the metadata file, you have to specify which column to be used for augmentation (i.e. corresponding column label). 


The outputs of the code are (written on the output directory)

1. Augmented data (Qiime2 FeatureTable[Frequency] artifact)
2. A copy of original data in the biom format (Qiime2 FeatureTable[Frequency] artifact)
3. Only for balancing: Metadata for augmented data (SampleData[ClassifierPredictions])
4. Only for balancing: A copy of the metadata for original data (SampleData[ClassifierPredictions])
5. The log file to keep track of progress and errors (optional, please look at parameters)


## Using TADA for augmentation
The training data size has a tremendous effect on the machine learning method performance. Generating new samples from the training data can be helpful. Using `feature-table.qza` and `phylogeny.qza` you can augment data (using hierarchy of binomials) to it using the following command

```
mkdir outputs_binom
qiime feature-engineering tada --i-phylogeny phylogeny.qza \
--i-otu-table feature-table.qza \
--o-orig-biom outputs_binom/original-feature-table.qza \
--o-augmented-biom outputs_binom/augmented-feature-table.qza \
--p-output-log-fp outputs_binom/logfile.log
```

Using the above code, the output files 

* `outputs_binom/logfile.log`
* `outputs_binom/augmented-feature-table.qza`
* `outputs_binom/original-feature-table.qza` 

should match 

* `binom/augmented-feature-table.qzaa`
* `binom/original-feature-table.qza` 

respectively. If you wish to use beta binomial generative model to generate new samples, you can use the following command instead

```
mkdir outputs_beta_binom
qiime feature-engineering tada --i-phylogeny phylogeny.qza \
--i-otu-table feature-table.qza \
--o-orig-biom outputs_beta_binom/original-feature-table.qza \
--o-augmented-biom outputs_beta_binom/augmented-feature-table.qza \
--p-output-log-fp outputs_beta_binom/logfile.log \
--p-stat-method beta_binom
```

Files 

* `outputs_beta_binom/augmented-feature-table.qza` 
* `outputs_beta_binom/original-feature-table.qza` 

should match 

* `beta_binom/augmented-feature-table.qza`
* `beta_binom/original-feature-table.qza` 

respectively.

### Using TADA for balancing datasets
In microbiome samples, the distribution of class labels (or cluster labels for unsupervised learning) is often unbalanced. This can cause overfitting and poor generalization of the machine learning method on new samples. You can use TADA to generate new samples for the underrepresented classes to make classes the same size. Using our example files `phylogeny.qza`, `feature-table.qza`, `metadata.csv`, and using labels in column `label` you can generate new samples to balance out datasets using 

```
mkdir outputs_binom_balancing
qiime feature-engineering tada --i-phylogeny phylogeny.qza \
--i-otu-table feature-table.qza \
--m-meta-data-file metadata.csv \
--m-meta-data-column label \
--o-orig-biom outputs_binom_balancing/original-feature-table.qza \
--o-augmented-biom outputs_binom_balancing/augmented-feature-table.qza \
--p-original-meta outputs_binom_balancing/original-metadata.csv \
--p-augmented-meta outputs_binom_balancing/augmented-metadata.csv \
--p-output-log-fp outputs_binom_balancing/logfile.log
```


This will generate the folder `outputs_binom_balancing`, and it will create the following files under this directory:

* `augmented-metadata.csv`: class/cluster labels of the new samples, first column: sample IDs, second column: labels
* `augmented-feature-table.qza`: generated features
* `logfile.log`: log file for generating features
* `original-feature-table.qza`: original features
* `original-metadata.csv`: original meta data file

The generated files 

* `outputs_binom_balancing/augmented-feature-table.qza`
* `outputs_binom_balancing/original-feature-table.qza`
* `outputs_binom_balancing/augmented-metadata.csv`
* `outputs_binom_balancing/original-metadata.csv` 

should match files available at 

* `binom_balancing/augmented-feature-table.qza`
* `binom_balancing/original-feature-table.qza`
* `binom_balancing/augmented-metadata.csv`
* `binom_balancing/original-metadata.csv` 

respectively. 

If you wish to use the Beta-Binomial generative model, you can use the following commands respectively

```
mkdir outputs_beta_binom_balancing
qiime feature-engineering tada --i-phylogeny phylogeny.qza \
--i-otu-table feature-table.qza \
--m-meta-data-file metadata.csv \
--m-meta-data-column label \
--o-orig-biom outputs_beta_binom_balancing/original-feature-table.qza \
--o-augmented-biom outputs_beta_binom_balancing/augmented-feature-table.qza \
--p-original-meta outputs_beta_binom_balancing/original-metadata.csv \
--p-augmented-meta outputs_beta_binom_balancing/augmented-metadata.csv \
--p-output-log-fp outputs_beta_binom_balancing/logfile.log \
 --p-stat-method beta_binom
```


The outputs are similar to what described above. Files 

* `outputs_beta_binom_balancing/augmented-feature-table.qza`
* `outputs_beta_binom_balancing/original-feature-table.qza`
* `outputs_beta_binom_balancing/augmented-metadata.csv`
* `outputs_beta_binom_balancing/original-metadata.csv` 

should match files available at 

* `beta_binom_balancing/augmented-feature-table.qza`
* `beta_binom_balancing/original-feature-table.qza`
* `beta_binom_balancing/augmented-metadata.csv` 
* `beta_binom_balancing/original-metadata.csv` 

respectively. 


The above command will generate enough number of samples for the least size cluster/class (in provided example, from group `1`) so that both groups have the same size. In this implementation of TADA, user can choose to continue generating samples so that the final size of each group be a multiple of initial size of the most frequent group. For example, if the most frequent gorup has `20` samples, and the least size group has `10` samples, and user wishes to have augmentation level of `5x`, then the final size of both classes will be `120 = (5)*20 + 20`. The augmentation level of `0x` means the user wants both classes the same size and no further augmentation (default). For example, the following command will peform a `5x` augmentation.

```
mkdir outputs_beta_binom_balancing_5x
qiime feature-engineering tada --i-phylogeny phylogeny.qza \
--i-otu-table feature-table.qza \
--m-meta-data-file metadata.csv \
--m-meta-data-column label \
--o-orig-biom outputs_beta_binom_balancing_5x/original-feature-table.qza \
--o-augmented-biom outputs_beta_binom_balancing_5x/augmented-feature-table.qza \
--p-original-meta outputs_beta_binom_balancing_5x/original-metadata.csv \
--p-augmented-meta outputs_beta_binom_balancing_5x/augmented-metadata.csv \
--p-output-log-fp outputs_beta_binom_balancing_5x/logfile.log \
--p-stat-method beta_binom \
--p-xgen 5
```

The outpus are similar to what described above. Please note that in the augmented meta data file, there are `110` samples with class `1` and `100` samples with class 0. Overal, you will have `120` samples for both classes. Also, the output files
 
* `outputs_beta_binom_balancing_5x/augmented-feature-table.qza`
* `outputs_beta_binom_balancing_5x/original-feature-table.qza`
* `outputs_beta_binom_balancing_5x/augmented-metadata.csv`
* `outputs_beta_binom_balancing_5x/original-metadata.csv`
 
should match files available at 

* `beta_binom_balancing_5x/augmented-feature-table.qza`
* `beta_binom_balancing_5x/original-feature-table.qza`
* `beta_binom_balancing_5x/augmented-metadata.csv`
* `beta_binom_balancing_5x/original-metadata.csv` 

respectively. 

Please note that, if you wish to have normalized counts, such thatcounts for each sample add up to one, use the flag `--p-normalized`. The default paramerter set outputs unnormalized counts.

## Seed number

Please note that you can specify the seed number for the random number generator  by passing `--seed [a number]`. By default, the seed number is `0`. 


## Import feature data into the 
You can import preprocessed feature table into Qiime using the following command

```
qiime tools import \
--input-path [path to feature table assuming BIOM v2.1.0] \
--output-path [ARTIFACT file path] \
--type FeatureTable[Frequency] \
--input-format BIOMV210Format
```

## Import phylogeny into the Qiime artifact
You can import the phylogeny stored as a __"newick"__ format into the Qiime artifact using the following command

```
qiime tools import \
--input-path [path to phylogeny in newick format]\
--output-path [ARTIFACT file path]\
--type Phylogeny[Rooted]
```

## Citation

* Erfan Sayyari, Ban Kawas, Siavash Mirarab, TADA: phylogenetic augmentation of microbiome samples enhances phenotype classification, Bioinformatics, Volume 35, Issue 14, July 2019, Pages i31–i40, [DOI](https://doi.org/10.1093/bioinformatics/btz394)
* The Treatment-Naive Microbiome in New-Onset Crohn’s Disease American Gut: an Open Platform for Citizen Science Microbiome Research Daniel McDonald, Embriette Hyde, Justine W. Debelius, James T. Morton, Antonio Gonzalez, Gail Ackermann, Alexander A. Aksenov, Bahar Behsaz, Caitriona Brennan, Yingfeng Chen, Lindsay DeRight Goldasich, Pieter C. Dorrestein, Robert R. Dunn, Ashkaan K. Fahimipour, James Gaffney, Jack A. Gilbert, Grant Gogul, Jessica L. Green, Philip Hugenholtz, Greg Humphrey, Curtis Huttenhower, Matthew A. Jackson, Stefan Janssen, Dilip V. Jeste, Lingjing Jiang, Scott T. Kelley, Dan Knights, Tomasz Kosciolek, Joshua Ladau, Jeff Leach, Clarisse Marotz, Dmitry Meleshko, Alexey V. Melnik, Jessica L. Metcalf, Hosein Mohimani, Emmanuel Montassier, Jose Navas-Molina, Tanya T. Nguyen, Shyamal Peddada, Pavel Pevzner, Katherine S. Pollard, Gholamali Rahnavard, Adam Robbins-Pianka, Naseer Sangwan, Joshua Shorenstein, Larry Smarr, Se Jin Song, Timothy Spector, Austin D. Swafford, Varykina G. Thackray, Luke R. Thompson, Anupriya Tripathi, Yoshiki Vázquez-Baeza, Alison Vrbanac, Paul Wischmeyer, Elaine Wolfe, Qiyun Zhu, The American Gut Consortium, Rob Knight
mSystems May 2018, 3 (3) e00031-18; [DOI](10.1128/mSystems.00031-18)
* Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. [DOI](https://doi.org/10.1038/s41587-019-0209-9)
