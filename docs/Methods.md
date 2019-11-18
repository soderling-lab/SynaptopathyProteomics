# Methods

## Mice
Shank2, Shank3, and Ube3a mice were gifts of Yong-Hui Jiang (Duke University).
Syngap1 were recieved from Gavin Rumbaugh (Scripps).

## Synaptosome Isolation.  
Synaptomes were prepared by sucrose gradient ultracentrifugation in a manner
derived from that of Gray and Whittaker (1962; see also
DiGiovanni et al., 2012 and Distler et al., 2014). 

Adult mice (aged 52-63 days old) were anesthetized with isoflurane and
euthanized by cervical decapitation. The brain was quickly removed and sliced
into 1mm sections (Zivic brain slicer; BSMAS001-1). Neocortex was removed and
flash frozen in liquid nitrogen. Striatum tissue was disected from slices
beween from -3.5 mm and -6.0 mm posterior to the anterior olfactory bulb and 
flash frozen. Tissue was stored in liquid nitrogen until 
synaptosome purification. 

Approximately 100-150 mg of tissue was homogenizing in 1mL isotonic sucrose
buffer (320 mM sucrose, 5 mM HEPES, 1mM EGTA, pH 7.4) containing protease
inhibitors (AEBSF 100mM, Sodium Orthovanadate 100mM, Leupeptin/peptstatin 2mg/mL). 
All steps were performed on ice or at 4*C and all buffers contained protease
inhibitors.
The crude homogenant was transfered into a 1.5 mL tube and centrifuged
for 10 min @ 1,000 g 4*C.  The supernatant containing the crude cytosolic fraction
was removed by pipette and transfered into a new 1.5 mL tube and centrifuged again 
for 20 min @ 12,000g 4*C. The supernatant was removed and then the crude synaptosome
pellet was resuspended in 500 uL resuspension buffer (320 mM sucrose, 
5mM Tris/HCl, pH 8.1). 

The resuspened crude synaptosome fraction was layered over a 0.8M-1.0M-1.2M 
discontinuos sucrose gradient, prepared from HEPES buffered sucrose in a 
Beckman Ultra-Clear 14x89mm tube (344059). 

Theis was centrifuged for 2 hours at 85,000 x g in SW 41 Ti Rotor.
Ultracentrifugation yielded a floating myelin fraction, a light membrane mixture
(0.8 M/1.0 M interface), a purified synaptosome fraction (1.0 M/1.2 M
interface), and a mitochondrial pellet. The purified synaptosome fraction was
removed by pipette, and then layer over a 10 mL 0.8 M sucrose cushion in a in a
fresh Ultra-Clear tube. This was spun 1 hour at 53,000 x g in SW 41 Ti Rotor to
pellet the purified synaptosomes. The supernatant was discarded and the
synaptosome pellete was stored at -80*C prior to processing for mass
spectrometry. 

## Sample Preparation for LC/MS. 
Get details from Greg.

## TMT Proteomics Raw peptide intensity data were exported from Proteome
Raw peptide data were exported from Proteome Discoverer (PD version 2.2, 
Thermo Fisher Scientific).
All data processing and analysis were performed in R and Python using custom written
scripts (see data and software availability). 

## Peptide level processing

### Sample loading normalization 
Sample loading normalization was performed to equalize total TMT-channel 
intensities within experimental batches. Sample loading normalization scales the
run-level intensity sums across all 11 channels within an experiment to be
equal. Intensity data were multiplied by a factor such that the mean sum of all columns
within an experiment were equal. This normalization step reflects the assumption
that an equal amount of protein was processed for MS from each TMT-replicate. 


### Peptide level filtering
Peptides that contained any missing QC values, more than two missing biological replicates, 
or more than 50% missing values overall within an experiment overall were removed.

### Peptide level imputing 
A small number of missing values were present at the peptide level 
(less than 1% of all values from both cortex and striatum). 
Missing values were inferred to be missing not at random (MNAR) due to their
left-shifted distribution compared to peptides without missing values.
Missing values were replaced using the k-nearest neighbors (KNN) 
algorithm as implemented by the impute.knn function from the impute package.

### QC filtering 
Each 11-plex TMT experiment included a common quality control (QC)
sample analyzed in triplicate. These QC samples were used to assess
intra-experimental variability as well as normalize between experiments.

Peptides that exhibit highly variable QC measurements were removed in a 
manner that was previously decribed. [Ping et al.,
2018](https://www.nature.com/articles/sdata201836).

Briefly, the three QC technical replicates within a batch
should yield identical measurements, and therefore the distribution of log2
ratios of the three QC channels for all peptides measured should fit a normal
distribution centered around zero. The standard deviation of this distribution
reflects technical variation in a peptide’s measurement and allows one to
evaluate a measurements precision. However, since variance differed based on
peptide mean intensity, we subdivided all QC ratios into 5 quantile bins. 
Peptides considered an outlier and removed if their mean QC
ratio exceeded a threshold of 4 standard devaiates from the mean of a bin. 

To identify peptides with highly variable QC measurements, peptide level QC data
were binned by intensity and peptides whose mean ratio (QC1/QC2, QC2/QC3,
QC1/QC3) were more than 4 standard deviations away from the mean the intensity
bin were considered outliers and removed.


## Protein level processing

### Protein summarization 
Proteins were summarized as the sum of all peptides mapped to a protein. 
Sample loading normalization was then performed across batches. 

### Intrabatch regression 
Each experimental cohort was prepared in two batches.
This was necessary because the ultra-centrifuge used held a maximum of six
samples. This intra-batch batch effect was recorded for 6/8 experiments.

This batch effect was quantified as the absolute value of the median biweight 
correlation between batch and a sample’s first principal component. 
A quantifiable batch effect was identified and removed
for Shank3, Syngap1, and Ube3a (cortex) and Shank2 and Ube3a (striatum)
experiments. If there was no evidence of a batch effect (bicor less than 0.1), then ComBat
was not applied. We corrected for this batch effect using the ComBat function from the sva package 
(Johnson et al., 2007). in two batches.


### IRS Normalization 
Internal reference scaling equalizes the protein-wise mean QC measurements 
across experimental batches.
We observed that although 82.7-87% (Cortex and Striatum, respectively) of
proteins were reliably quantified across all four experiments, peptide
identification overlap was considerably lower. Only 24-34% of peptides were
quantified across all four experiments (striatum, cortex). That is, although
~80% of proteins are identified across all four experiments, these proteins are
identified by different peptides. This is due to random sampling of peptides at
the MS2 level. To account for this bias, we utilized our QC samples to perform
internal reference standard (IRS) normalization. IRS normalization scales QC
measurements to be equal, normalizing protein measurements between batches
(Plubell et al., 2017).

As IRS normalization assumes that IRS samples
are identical, outlier QC samples were identified iteratively and removed using
the method of Oldham et al., (2008). From cortex, two QC outliers (threshold =
Z-score Normalized Connectivity -2.5) were removed: Abundance: F2: 129N, Sample,
QC, 38824, Ube3a; Abundance: F1: 129N, Sample, QC, 36479, Syngap1 QC. No QC
outliers were identified in the striatum experiments. 


### Protein level filtering and imputing 
Proteins that were identified by a single peptide or identified 
in less than 50% of all samples were removed (Cortex = 94, 15; Striatum = 119,10). 
Any remaining missing values were inferred to be MNAR and imputed with the KNN
algorithm.


### Combining Cortex and Striatum data 
Data from cortex and striatum were combined, keeping proteins identified in both tissues. 
We developed a normalization approach, TAMPOR normalization, to combine data from this tissues. TAMPOR
normalizes data relative to WT samples.

A Tunable Approach for Median-Polish of Ratio (TAMPOR). Briefly, this
normalization method scales WT samples from cortex and striatum to be equal.
Only proteins that were identified in both cortex and striatum datasets were
used (N=3,022). Within a batch, biological replicates are multiplied by the same
scaling factor, so this approach does not affect intra-batch comparisons.
Finally, prior to differential expression and WGCNA analysis, samples with low
sample z-score normalized network connectivity (Z.Ki less than −2.5) (outlier samples)
were identified iteratively and removed Oldham et al., 2008. One cortex sample
was identified as an outlier and removed. 

### Sample level outliers Outlier samples were identified using the method of
Oldham et al.,.  I sample outlier was identified (Syngap1 HET b7.131N).

In all, 2,918 proteins were reliably quantified across 63 biological samples.

## Differentially abundant proteins We utilized a generalized linear model to
compare protein abundance in pooled tissue specific WT and KO samples as
implemented by the glmQLFit _edgeR_ package (Robinson et al., 2010). 
A general linear model was fit to samples. glmQLFit function.  Differences
between groups were evaluated using the glmQLFTest function.  P value correction
Benjamini Hochberg.  Significant if p.value is less than 0.05.

## Coexpression network construction A WT and KO coexpression network were
constructed using the bicor function from the WGCNA function. Bicor is a robust
alternative to Spearman correlation.

Define an adjacency matrix as the mxm correlation matrix of m genes.  A
correlation between gene m1 and mi is a measure of the distance between the two
genes. The bicor correlation metric was used as a robust alternative to the
Pearson correlation.

The quality of a network partition can be described in terms of modularity,
roughly how dense clusters are and well seperated...  Leiden algorithm performs
clustering with optimization of a quality function like modularity.

## Clustering

### WGCNA WGCNA is a heirarchical clustering approach used to identify modules
or clusters of co-expressed genes or proteins.

The WGCNA function is expensive to evaluate and depends upon a large number of
parameters. We optimized the parameters of the WGCNA function using Beyesian
hyperparameter optimization using the gp_minimize function from the
[Skikit-Optimize](https://scikit-optimize.github.io/) Python module.

The optimization of modularity cannot be performed by an exhaustive search as
the number of different partitions of the graph grows exponentially with the
number of nodes N (Brandes, Delling... 2008).

Quality of the co-expression graph was assesed using the modularity of signed
weighted networks (REF).  Minimize the inverse of Qws. 

Optimizing these parameters resulted in a network with 3 large clusters. This
result is similar to other modularity-optimization approaches.

This is result is consistent with the resolution limit, an acknowledge
limitation of modularity optimization based clustering approaches.

To overcome this limitation, several different groups have added a resolution
parameter...

### Leiden algorithm The Leiden algorithm was used to perform meso-resolution
clustering of the co-expression networks at 100 resolutions from 0-1.  As
implemented by the Python leidenalg package.

We focused on partitions with the most biological information as assesed by
summing the log10(pvalues) for enriched processes for every module.

Thus we identified 1 optimized WT and 1 optimized KO resolution.

## Module Sel-preservation Module quality was ensured by permutation testing.
Modules with any insignificant quality statistics were removed.

## Comparing WT and KO networks
In order to identify modules that are divergent in WT and KO coexpression graphs,
we again utilized a permutation approach. (Ritchie et al., 2016).

A permutation procedure was employed to characterize the distribution of
each statistical test under the null hypothesis of non-replication and non-preservation. Specifically, each module preservation statistic was re-calculated
when shuffling the node labels in the test dataset. The node labels in the discovery dataset were left unchanged. Nodes that were not present in both the
discovery and test dataset were ignored both when calculating the module
preservation statistics and when shuffling the node labels in the test dataset.
Under the alternate hypothesis of replication/preservation, the test statistics
calculated on the non-permuted dataset were expected to be higher than
when calculated on random sub-graphs in the test dataset. Permutation
p values were then calculated from these null distributions using the estimator
described by Phipson and Smyth (2010), which provides a conservative estimate of the p value appropriate for multiple testing adjustment (Supplemental
Experimental Procedures).

## Data and Software availability
Data processing and statistical analyses were performed using R version 3.5.1
with custom written scripts. Plots were generated with the ggplot2 package
(Wickham 2016). A complete list of dependencies is given in supplemental file X.
All data and code are available online
(https://github.com/twesleyb/SynaptopathyProteomics). This private repository
will be made public upon publishing. 

# Acknowledgments
We are indebted to Duke’s proteomics facility staff for their efforts in the
data collection and preliminary data analysis. We would like to thank Dr. Gavin
Rumbaugh for his generous donation of the Syngap1 mice used in this study. We
would also like to thank Dr. Peter Mucha for helpful discussions about the
analysis. 


## module comparison notes.

# For example:
# Given a network partition at a given resolution:
# Get observed statistics for a module in the opposite (test) graph (e.g. KO module avg.edge weight in WT graph).
# Get nulls from 10,000 randomizations of the test graph.
# Compare observed versus NULL distributions.
# Correct p.values for n comparisons.
# If observed statistic is significantly greater than null -> preserved.
# If observed statistic is significantly less than null --> divergent.
# WT Modules that are divergent in the KO graph coorespond to LOF.
# KO Modules that are divergent in the WT graph coorespond to GOF.

# Criterion for preservation/divergence can be weak or strong.
# If weak, then requirement is ANY significant stats.
# If strong, then requirement is ALL significant stats.
# User can select which stats to use.

# Possible permutations:
# 1. enforce self-preservation
#    * strong preservation -- all stats - NO DIVERGENT WT MODULES
#    * weak preservation -- all stats
#    * strong preservation -- 2 stats - NO DIVERGENT WT MODULES?

# 2. no self-pres --> remove modules smaller than 5 proteins
#    * strong preservation -- all stats
#    * weak preservation -- all stats
#    * strong preservation -- 2 stats *** THIS seems like the prefered option ***

# 3. Single resolution : best biological "signal" - self-pres with TOM -> recalculate best resolution.
#    * strong preservation -- all stats
#    * weak preservation -- all stats
#    * strong preservation -- 2 stats


# Module statistics:
# 1. avg.weight - (average edge weight) assumes positive edges - calculated 
#    from network.
# 2. coherence (module coherence) - calculated from modules summary profile.
# 3. cor.cor (concordance of correlation structure) - calculated from
#    correlation matrix.
# 4. cor.degree - assumes positive edges - calculated from network.
# 5. cor.contrib (concordance of node contribution) - calculated from modules 
#    summary profile. 
# 6. avg.cor (density of correlation structure) - calculate from correlation 
#    matrix.
# 7. avg.contrib (average node contribution) - calculated from modules summary
