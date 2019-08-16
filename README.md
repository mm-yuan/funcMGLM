
# funcMGLM: fit Multivariate Generalized Linear Models with Functional Forms

Provides functions that (1) fit classical multivariate GLM, which the correlations among multiple dependent variables are expressed in the covariance matrix
(2) fit multivariate GLM with functional forms(underlying value, slope, AUC), which provides direct interpretation for the correlations

## Usage

Based on the PBC dataset:

```
library(JMbayes)
library(MASS)
library(splines)
library(matrixStats)

fit = mvglmer(formulas = list(log(serBilir) ~  ns(year,2) + sex + drug + (ns(year,2)| id), 
                              spiders ~ year + albumin + sex + (year| id)),
             data = pbc2, families = list(gaussian, binomial),
             optionHC = "HC", scaling = "standardize")
```


## Installation


