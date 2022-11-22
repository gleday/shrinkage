# shrinkage
An R package to fit Bayesian regression models with shrinkage priors.

Linear regression models with global shrinkage priors:
* normal-[inverse-Gamma](https://en.wikipedia.org/wiki/Inverse-gamma_distribution)
* normal-[beta prime](https://en.wikipedia.org/wiki/Beta_prime_distribution)
* normal-[inverse-Gaussian](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution)
* normal-[Gamma](https://en.wikipedia.org/wiki/Gamma_distribution)
* normal with [empirical Bayes](https://en.wikipedia.org/wiki/Empirical_Bayes_method) (closed-form inference)

Linear regression models with local shrinkage priors:
* normal-Gamma ([Griffin and Brown, 2010](https://projecteuclid.org/euclid.ba/1340369797))
* normal-inverse-Gaussian ([Caron and Doucet, 2008](http://doi.acm.org/10.1145/1390156.1390168)) 
* normal-Beta-Prime ([Bai and Gosh, 2019](http://www3.stat.sinica.edu.tw/ss_newpaper/SS-2019-0037_na.pdf))
* accommodate shrinkage for grouped variables

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

