# SignatureLikelihoods
Mutational signature analysis for low statistics SNV data

Functionalities:
* Construct SNV matrix with trinucleotide context
* Calculates cosine similarity and likelihood
* Generate MC based on single signature probability distributions as well as 
* Tunes the cutoff on these parameters to maximize sensitivity with a reasonable false positive rate
* Plotting macros
  * Heat map for cosine similary and likelihood
  * Plots the change in sensitivity and false positive rate by changing the threshold on sensitivity/cosine similarity
  * Plots the sensitivity and false positive rate as a function of number of SNVs for tuned values of cutoff

An example on how to use the package can be found at example.R

