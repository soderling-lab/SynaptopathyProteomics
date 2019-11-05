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
Sample loading normalization was performed to equalize run level intensity 
sums within experimental batches.

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
by inspection of their density plot and imputed using the KNN algorithm as 
implemented by the `impute.knn` funciton from the _impute_ package.

### QC filtering
Each 11-plex TMT experiment included a quality control (QC) sample 
analyzed in triplicate. These QC samples can be used to asses 
intra-experimental variability.

Peptides that exhibit highly variable QC measurements were filtered.

To identify peptides with highly variable QC measurements,
the peptide level QC data were binned by intensity and peptides 
whose mean ratio (QC1/QC2, QC2/QC3, QC1/QC3) were more than 4 standard
deviations away from the mean the intensity bin were considered outliers 
and removed.

Filter peptides based on QC:
Cortex:
Striatum: 182, 67, 73, 75

## Protein level processing

### Protein summarization
Summarize to protein level by summing peptides for all proteins.

## Intrabatch regression (ComBat)
Cortex: Shank3, Syngap1,Ube3a
Striatum: Shank2, Ube3a 

Each experimental cohort was prepared in two batches.
This was necessary because the ultra-centrifuge used held a maximum of
siz samples. This intra-batch batch effect was recorded for 6/8 experiments.
The `ComBat` function from the _sva_ package was used to remove this intra-batch
batch effect before correcting for inter-batch batch effects between batches
with IRS normalization. Evidence of a batch effect was taken to be cor(PCA,batch)>0.1.

## IRS Normalization
Internal reference scaling equalizes the protein-wise means of reference (QC)
samples across all batchs. the IRS normalization accounts for the random 
sampling of peptides at the MS2 level whic results in the identifiaciont of
proteins by different peptides in each experiemnt. (Plubell et al., 2017).

## QC outlier removal
Cortex: 2
Striatum: 0
IRS normalization utilizes QC samples as reference samples, outlier QC 
measurements caused by interfeerence or other technical artifacts would
biase this step. QC sample outliers were identified in the maner of Oldham
et al., 2016 and removed if Z-score normalized sample connectivity was
less than 2.5.

## Protein level filtering and imputing.
Proteins identified by a single peptide as well as proteins identified in
less than 50% of samples removed.
Any remaining missing values were imputed using th KNN algorithm.

Cortex: 52 protein identified by single peptide, 1 less than 50% sampels.
Striatum: 80 prteins identified by singple peptide, 4 proteins in less than 50%.

## Combining Cortex and Striatum data
TAMPOR normalization

