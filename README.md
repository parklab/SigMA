# SignatureLikelihoods
SigMA is a signature analysis tool optimized to detect the mutational signature associated to HR defect, Signature 3, from hybrid capture panels, exomes and whole genome sequencing. For panels with low SNV counts, conventional signature analysis tools do not perform well while the novel approach of SigMA allows it to detect Signature 3-positive tumors with 74% sensitivity at 10% false positive rate.

One novelty of SigMA is a likelihood based matching: We associate a new patient's mutational spectrum to subtypes of tumors according to their signature composition. The subtypes of tumors are defined using the WGS data from ICGC and TCGA consortia, by a clustering of signature fractions with hierarchical clustering. The likelihood of the sample to belong to each tumor subtype is calculated, and the likelihood of Signature 3 is the sum of the likelihoods of all Signature 3-positive tumor subtypes. The second novel step is the multivariate analysis with gradient boosting machines, which allows us to obtain a final score for presence of Signature-3 combining likelihood with cosine similarity and exposure of Signature 3 obtained with non-negativel-least-squares (NNLS) algorithm. The multivariate analysis allows us to automatically handle different sequencing platforms. For different platforms different methods for signature analysis become more efficient, e.g. for WGS data it is not necessary to associate the tumor to a subtype of tumors, because it is possible to determine Signature 3 with NNLS acurately. We have a new feature also for these cases and we calculate the likelihood of NNLS decomposition to be unique. This likelihood value was found to be the most influential feature in the multivariate analysis.

### Functionalities:

* Construct SNV matrix with trinucleotide context
* Calculates cosine similarity and likelihood
* Generate MC based on single signature probability distributions as well as 
* Tunes the cutoff on these parameters to maximize sensitivity with a reasonable false positive rate
* Plotting macros
  * Heat map for cosine similary and likelihood
  * Plots the change in sensitivity and false positive rate by changing the threshold on sensitivity/cosine similarity
  * Plots the sensitivity and false positive rate as a function of number of SNVs for tuned values of cutoff

An example on how to use the package can be found at example.R

