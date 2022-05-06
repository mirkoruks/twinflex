# twinflex
twinflex is a wrapper function for OpenMx. It allows you to estimate complex twin models with just one line of code. Furthermore you can combine different types of twin models (univariate, GxE, with or without covariates, Cholesky vs. ACE-beta) individually. It is still a beta version. So if you find any bugs or have any questions, please send me a message (mirko.ruks@uni-bielefeld.de) or start an issue here on Github. 

# Installation
In order to install the package from GitHub use:

```
devtools::install_github('mirkoruks/twinflex')
```


# Documentation
In the next time I will add more and more to the documentation, including coding examples.

# Output Interpretation
At the moment, the following types of parameters are estimated: 1) path coefficients, 2) variances and covariances, 3) means and intercepts, 4) thresholds

## Path coefficients
Here, the notation is simple. The name of a path coefficient consists of two pars: The label and the suffix. The label is a string. The suffix is numeric. The label depends on the type of the effect: 

- Additive genetic effects have the label "a"
- Non-additive genetic effects have the label "d"
- Shared environmental effects have the label "c"
- Non-shared environmental effects have the label "e"
- Phenotypic effects ("beta" effects) have the label "b"
- Random Long-formatted covariates (part of the covariance matrix) have the label "bcovl"
- Random Wide-formatted covariates (part of the covariance vector) have the label "bcovw"
- Fixed Wide- and Long-formatted covariates (part of the mean vector) have the label "bcov" 

The suffixes consists of two numerics. The first refers to the position of the dependent variable in the vector given in the `acevars` argument. The second refers to the position of the independent variable in the vector given in the `acevars` or `covvars` argument. So, e.g. a21 refers to the effect of the genetic component of acevar no. 1 on acevar no. 2.

Interaction effects are somewhat longer. There is a label for all interaction effects "bm", then follows an integer referring to the number of moderator (at the moment it is possible to use 5 different moderators. However, that you are able to use 5 does not imply that you should use 5.). Then comes the "normal" name of the path which is moderated. For example bm1a21 refers to the interaction of the effect of the genetic component of acevar no. 1 on acevar no. 2 by moderator no. 1.

## Variances and Covariances
## Means and Intercepts
## Thresholds

## Outlook
I'm sure I will add a feature that shows the output "in words". For example, for a model with `acevars = c('IQ','yeduc')`, instead of a21 it will give something like A(IQ) -> yeduc. A phenotypic effect would just be X -> Y. A moderation would be something like X -> Y BY M
