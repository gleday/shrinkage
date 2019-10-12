# shrinkage
An R package implementing various regression models with shrinkage priors.

Regression models with global shrinkage priors:
* Gaussian (ridge) prior with global variance endowed with [inverse-gamma](https://en.wikipedia.org/wiki/Inverse-gamma_distribution), [gamma](https://en.wikipedia.org/wiki/Gamma_distribution), [beta prime](https://en.wikipedia.org/wiki/Beta_prime_distribution) or [inverse-Gaussian](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution) or estimated using an [empirical Bayes method](https://en.wikipedia.org/wiki/Empirical_Bayes_method)
* Gausian prior for grouped variables (group ridge)

Regression models with local shrinkage priors:
* normal-Gamma ([Griffin and Brown, 2010](https://projecteuclid.org/euclid.ba/1340369797))
* normal-inverse-Gaussian ([Caron and Doucet, 2008](http://doi.acm.org/10.1145/1390156.1390168)) 
* normal-Beta-Prime ([Bai and Gosh, 2019](http://www3.stat.sinica.edu.tw/ss_newpaper/SS-2019-0037_na.pdf))

Regression models with global-local shrinkage priors:
* horseshoe prior (hs)
* R2-D2 prior (r2d2)
* Dirichlet-Laplace prior (dl)

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

