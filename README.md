# SigMoS
Model selection for robust learning of mutational signatures

Our package provides functions for estimating mutational signatures from mutational count data from cancer patients. 
Our analysis consists of the following steps:

1. apply Poisson-NMF (cross-validation for choosing number of signatures assuming the Poisson distribution)

2. Goodness of fit: residual analysis to validate Poisson assumption.

3. If overdispersion is present: NB-NMF (cross-validation for choosing number of signatures assuming the Negative Binomial distribution) and its patient specific counterpart NB$_\text{N}$-NMF (dispersion parameter is patient specific)

4. Likelihood ratio test to choose between NB-NMF and NB_N-NMF

5. Residual analysis of the final Negative Binomial NMF model

See our vignettes for more details.

Our preprint can be found at https://arxiv.org/abs/2206.03257

# Examples on how to use the functions 

The package can be installed using the package **devtools**

```{r}
devtools::install_github("MartaPelizzola/SigMoS")
library(SigMoS)
```

Assuming Poisson as the underlying distribution for NMF we can apply our model selection procedure:

```{r poisson}
k_vec <- 2:8                          # a range of the number of signatures assumed
nsig_pois <- list()                   # list of outputs for each assumed number of signatures
cost <- c()                           # vector of cost values from the model selection procedure
for (i in 1:length(k_vec)){
  nsig_pois[[i]] = sigmos(data = BRCA21, k = i, method = "Poisson")
  cost = c(cost, nsig_pois[[i]]$cost_k)
}
k_poisson = k_vec[which.min(cost)]    # the assumed number of signatures with the smallest cost
```

With the resulting optimal number of signatures we apply Poisson NMF to obtain the signature and exposure matrices under the Poisson distribution:

```{r poissonNMF}
res_pois = NMFPois(M = BRCA21, N = k_poisson)
```

Given these results we can check the residuals and obtain the following residual plot:

```{r poissonRSD, warning = FALSE}
W_p <- res_pois$E
H_p <- res_pois$P
rsd_p = (BRCA21 - H_p%*%W_p)
rsd_norm = rsd_p/(sqrt(H_p%*%W_p))
funupper = function(x) 2*sqrt(x)
funlower = function(x) -2*sqrt(x)
dataPois = data.frame(obs = as.vector(H_p%*%W_p),Residuals = as.vector(rsd_p), residual_norm = as.vector(rsd_norm))
ggplot(dataPois, aes(x = obs, y = Residuals))+
  geom_point()+
  scale_x_log10()+
  stat_function(fun = funupper, color = "#0072B2", size = 1 ) + 
  stat_function(fun = funlower, color = "#0072B2", size = 1 ) +
  ylim(-100,100) + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(size=14),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        line = element_line(margin(b = 6, l = 6)), 
        legend.position = "none", 
        plot.title = element_text(size = 20)) +
  xlab("Estimated mean") + ylab("Raw residuals")
```

Following the same procedure as for the Poisson case, we first find the optimal number of signatures under the Negative Binomial model with our model selection approach and then we estimate the signature and exposure matrices with NMF.


```{r NBpatientspecific, warning = FALSE}
k_vec <- 2:8
nsig_nb_ps <- list()
cost <- c()
for (i in k_vec){
  nsig_nb_ps[[i-1]] = sigmos(data = BRCA21, k = i, method = "NB", patient_specific = TRUE)
  cost <- c(cost, nsig_nb_ps[[i-1]]$cost_k)
}
k_NB_ps <- k_vec[which.min(cost)]
alpha_ps <- alphaNR(BRCA21,k_NB_ps, patient_specific = TRUE)
res_nb_ps = NMFNB(M = BRCA21, N = k_NB_ps, alpha = alpha_ps)
```
The residual plots for the negative binomial can be generated using the same code as the one for the Po-NMF.


# Authors
Ragnhild Laursen (ragnhild@math.au.dk)

Marta Pelizzola (marta@math.au.dk)

