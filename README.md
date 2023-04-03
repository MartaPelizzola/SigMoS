# SigMoS package 

The R package SigMoS provides functions to extract signatures and exposures from mutational count data from cancer patients. 

We provide functions to extract signatures and exposures based on the Negative Binomial and Poisson distribution and a function to perform model selection to estimate the number of signatures to be reconstructed. 
The package also contains a function to compare signatures based on their cosine similarity and a function to estimate the patient specific dispersion parameters for the Negative Binomial model.

For more information about SigMoS see our manuscript on bioRxiv "Model selection and robust inference of mutational signatures using Negative Binomial non-negative matrix factorization" by Pelizzola et al.  
=======

A typical mutational count data set we considered for the SigMoS package is an R matrix with 96 rows (equal to the number of mutation types with 6 base mutations, when assuming strand symmetry and the 4 flanking nucleotides on each side) and as many columns as many patients are included in the data set. 

In order to obtain a mutational matrix from a data file different functions can be used according to the format of the data file. Some examples are:

```{r readdata}
library(vcfR)
?read.vcfR
library(BiocManager)
BiocManager::install("maftools")
library(maftools)
?read.maf
```

We show below the different options available in this R package. 

## 1. SigMoS with Poisson-NMF
When choosing the Poisson distributional assumption we can apply the model selection function with the Poisson model as follows:
```{r poisson}
k_vec <- 2:8
nsig_pois <- list()
cost <- c()
for (i in k_vec){
  nsig_pois[[i-1]] = sigmos(data = BRCA21, k = i, method = "Poisson")
  cost = c(cost, nsig_pois[[i-1]]$cost_k)
}
k_poisson = k_vec[which.min(cost)]
```

## 2. SigMoS with NB-NMF
Following the same procedure as for the Poisson case, we can also choose the Negative Binomial distributional assumption and find the optimal number of signatures under the Negative Binomial model with our model selection approach.

```{r NBpatientspecific, warning = FALSE}
k_vec <- 2:8
nsig_nb_ps <- list()
cost <- c()
for (i in k_vec){
  nsig_nb_ps[[i-1]] = sigmos(data = BRCA21, k = i, method = "NB", patient_specific = TRUE)
  cost <- c(cost, nsig_nb_ps[[i-1]]$cost_k)
}
k_NB_ps <- k_vec[which.min(cost)]
```
See our vignettes for more details by typing:
```{r vignette, warning = FALSE}
browseVignettes("SigMoS") 
```
A working example of the different functions offered in our package can be found at: 

# Authors
Ragnhild Laursen (ragnhild@math.au.dk)

Marta Pelizzola (marta@math.au.dk)
