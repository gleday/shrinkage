# shrinkage
An R package to fit Bayesian regression models with shrinkage priors.

Linear regression models with global shrinkage priors:
* normal-[inverse-Gamma](https://en.wikipedia.org/wiki/Inverse-gamma_distribution)
* normal-[beta prime](https://en.wikipedia.org/wiki/Beta_prime_distribution) (which includes normal–half Cauchy)
* normal-[inverse-Gaussian](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution)
* normal-[Gamma](https://en.wikipedia.org/wiki/Gamma_distribution)
* normal with [empirical Bayes](https://en.wikipedia.org/wiki/Empirical_Bayes_method) (closed-form inference)

Linear regression models with local shrinkage priors, which accommodate grouped variables and some sparse-inducing priors:
* normal-inverse-Gamma ([Armagan and Zaretzki, 2010](https://link.springer.com/article/10.1007/s00180-010-0186-4))
* normal-Gamma ([Griffin and Brown, 2010](https://projecteuclid.org/euclid.ba/1340369797))
* normal-inverse-Gaussian ([Caron and Doucet, 2008](http://doi.acm.org/10.1145/1390156.1390168)) 
* normal-Beta-Prime ([Bai and Gosh, 2019](http://www3.stat.sinica.edu.tw/ss_newpaper/SS-2019-0037_na.pdf))

## Installation

To install **shrinkage** from R:

```R
# Install/load R package devtools
install.packages("devtools")
library(devtools)

# Install/load R package beam from github
install_github("gleday/shrinkage")
library(shrinkage)
```

## Vignettes

Vignette describing in details the implemented algorithms:
```R
vignette("Algorithms", package = "shrinkage")
```
Vignette providing code for fitting models using [Stan](https://mc-stan.org)
via the R package [rstan](https://cran.r-project.org/package=rstan):
```R
vignette("Stan", package = "shrinkage")
```
