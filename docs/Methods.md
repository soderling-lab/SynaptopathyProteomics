# TMT Data Preprocessing

## Raw data
Raw peptide intensity data were exported from Proteome Discoverer 
(PD version 2.2, Thermo Fisher Scientific).

From cortex, 43,651 unique peptides cooresponding to 3,305 unique proteins
were identified.

From striatum, 42,288 unique peptides cooresponding to 3,493 unique proteins 
were identified. 

## Peptide level processing

### Sample loading normalization
Sample loading normalization was performed to equalize total TMT-channel 
intensity within experimental batches.

Intensity data were multiplied by a factor such that the mean sum of all 
columns within an experiment were equal. This normalization step 
reflects the assumption that an equal amount of protein was processed for
MS from each TMT-replicate. 

### Peptide level filtering
Proteins with missing QC replicates, more than 2 missing biological
replicates, or more than 50% missing values within an experiment overall
were removed.

### Peptide level imputing
Any remaining missing values were infered to be missing at random
by inspection of the density plots of peptides with missing and non-
missing values and imputed using the KNN algorithm as 
implemented by the `impute.knn` function from the _impute_ package.

### QC filtering
Each 11-plex TMT experiment included a quality control (QC) sample 
analyzed in triplicate. These QC samples were used to assess
intra-experimental variability as well as normalize between experiments.

Peptides that exhibit highly variable QC measurements were removed in a
manner similar to that developed by [Ping et al., 2018](https://www.nature.com/articles/sdata201836).

QC replicates should yield identical measurements. Thus the ratio of
QC samples should be equal to 1. 

To identify peptides with highly variable QC measurements,
peptide level QC data were binned by intensity and peptides 
whose mean ratio (QC1/QC2, QC2/QC3, QC1/QC3) were more than 4 standard
deviations away from the mean the intensity bin were considered outliers 
and removed.

## Protein level processing

### Protein summarization
Proteins were summarized as the sum of all their peptides.

### Intrabatch regression
Each experimental cohort was prepared in two batches.
This was necessary because the ultra-centrifuge used held a maximum of
six samples. This intra-batch batch effect was recorded for 6/8 experiments.
The `ComBat` function from the _sva_ package was used to remove this intra-batch
batch effect before correcting for inter-batch batch effects between batches
with IRS normalization. Evidence of a batch effect was taken to be bicor(PCA,batch)>0.1.
A quantifiable batch effect was identified and removed for Shank3, Syngap1, and Ube3a
(cortex) and Shank2 and Ube3a (striatum) experiments.

### IRS Normalization
Internal reference scaling equalizes the protein-wise means of reference (QC)
samples across all batches. IRS normalization accounts for the random 
sampling of peptides at the MS2 level whic results in the identifiaciont of
proteins by different peptides in each experiemnt (Plubell et al., 2017).

### QC outlier removal
IRS normalization utilizes QC samples as reference samples, outlier QC 
samples caused by interference or other technical artifacts would
biase this step. QC sample outliers were identified in the maner of Oldham
et al., 2016 and removed if Z-score normalized sample connectivity to all
other QC samples was less than 2.5. Two QC outlier samples were identified in
cortex samples and removed.

### Protein level filtering and imputing
Proteins identified by a single peptide as well as proteins identified in
less than 50% of samples were removed.

Any remaining missing values were imputed using the KNN algorithm for missing
at random data.

### Combining Cortex and Striatum data
Finally we combined cortex and striatum data, keeping proteins identified in both
tissues. We developed a normalization approach, TAMPOR normalization, to combine
data from this tissues. TAMPOR normalizes data relative to WT samples.

In all, XXX proteins were reliably quantified across YYY biological samples.

## Differentially abundant proteins
We utilized a generalized linear model to compare protein abundance in
pooled tissue specific WT and KO samples as implemented by the glmQLFit
_edgeR_ package

A general linear model was fit to samples. glmQLFit function.
Differences between groups were evaluated using the glmQLFTest function.
P value correction Benjamini Hochberg.
Significant if p.value is less than 0.05.

## Coexpression network construction
A WT and KO coexpression network were constructed using the bicor 
function from the WGCNA function. Bicor is a robust alternative to 
Spearman correlation.

Define an adjacency matrix as the mxm correlation matrix of m genes.
A correlation between gene m1 and mi is a measure of the distance between
the two genes. The bicor correlation metric was used as a robust alternative
to the Pearson correlation.

The quality of a network partition can be described in terms of modularity,
roughly how dense clusters are and well seperated...
Leiden algorithm performs clustering with optimization of a quality function
like modularity.

## Clustering

### WGCNA
WGCNA is a heirarchical clustering approach used to identify
modules or clusters of co-expressed genes or proteins.

The WGCNA function is expensive to evaluate and depends upon
a large number of parameters. We optimized the parameters of the
WGCNA function using Beyesian hyperparameter optimization using 
the gp_minimize function from the Skikit-Optimize Python module.

The optimization of modularity cannot be performed by an exhaustive search
as the number of different partitions of the graph grows exponentially with the
number of nodes N (Brandes, Delling... 2008).

Quality of the co-expression graph was assesed using the modularity 
of signed weighted networks (REF).
Minimize the inverse of Qws. 

Optimizing these parameters resulted in a network with 3 large
clusters. This result is similar to other modularity-optimization 
approaches.

This is result is consistent with the resolution limit, an acknowledge
limitation of modularity optimization based clustering approaches.

To overcome this limitation, several different groups have 
added a resolution parameter...

### Leiden algorithm
The Leiden algorithm was used to perform meso-resolution clustering of
the co-expression networks at 100 resolutions from 0-1.
As implemented by the Python leidenalg package.

We focused on partitions with the most biological information as 
assesed by summing the log10(pvalues) for enriched processes for 
every module.

Thus we identified 1 optimized WT and 1 optimized KO resolution.

## Module Sel-preservation
Module quality was ensured by permutation testing. Modules with 
any insignificant quality statistics were removed.

## Comparing WT and KO networks

