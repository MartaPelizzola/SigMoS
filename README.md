# SigMoS
Cross-validation and model selection for robust learning of mutational signatures

Our analysis consists of the following steps:

1. apply Poisson-NMF (cross-validation for choosing number of signatures assuming the Poisson distribution)

2. Goodness of fit: residual analysis to validate Poisson assumption.

3. If overdispersion is present: NB-NMF (cross-validation for choosing number of signatures assuming the Negative Binomial distribution) and its patient specific counterpart NB$_\text{N}$-NMF (dispersion parameter is patient specific)

4. Likelihood ratio test to choose between NB-NMF and NB_N-NMF

5. Residual analysis of the final Negative Binomial NMF model

See our vignettes for more details

Our preprint can be found at 

# Authors
Ragnhild Laursen (ragnhild@math.au.dk)
Marta Pelizzola (marta@math.au.dk)

