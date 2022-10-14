<<<<<<< HEAD
# SigMoS package 

This package includes a range of different functions.
=======
# SigMoS
A typical mutational count data set we considered for the SigMoS package has 96 rows (equal to the number of mutation types with 6 base mutations, when assuming strand symmetry and the 4 flanking nucleotides on each side) and as many columns as many patients are included in the data set. 

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
See our vignettes for more details

# Authors
Ragnhild Laursen (ragnhild@math.au.dk)

Marta Pelizzola (marta@math.au.dk)

>>>>>>> d7805b65e1cea0e44490076ae5b7d74603995008
