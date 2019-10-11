# shrinkage
An R package implementing various regression models with shrinkage priors, including:
* Gaussian (ridge) prior with global variance endowed with [inverse-gamma](https://en.wikipedia.org/wiki/Inverse-gamma_distribution), [gamma](https://en.wikipedia.org/wiki/Gamma_distribution), [beta prime](https://en.wikipedia.org/wiki/Beta_prime_distribution) or [inverse-Gaussian](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution)
* group ridge (bgr)

To do's:
* horseshoe prior (hs)
* R2-D2 prior (r2d2)
* normal-Beta-Prime (nbp)
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

