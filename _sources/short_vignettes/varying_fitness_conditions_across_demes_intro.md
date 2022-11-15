---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(varying_fitness_conditions_across_demes_intro)=

# Introduction

The following vignettes describe how to set up models where the details of genetic value/fitness calculations differ between demes.

## Different effect sizes in different demes.

The next several vignettes describe how to set up scenarios where the effect size of a mutation differs between demes.
Briefly, simulating such models entails:

* Using {class}`fwdpy11.mvDES` to define a multivariate distribution of effect sizes.
* Passing new parameters to our genetic value types to indicate how many demes to expect during the simulation.

(correlation_and_covariance_matrices)=

## Correlation and covariance matrices

The next sections describe types that require covariance matrices as inputs.
Many people may find it easier to think about correlations between variables than covariances.

This section shows how to start with a *correlation* matrix and convert it to a *covariance* matrix.

Starting with a 3x3 correlation matrix:

```{code-cell} python
import numpy as np
corr = np.matrix([1,0.25,0.9,0.25,1,0.5,0.9,0.5,1]).reshape(3,3)
print(corr)
```

We wish to convert this to a covariance matrix with the following standard deviations:

```{code-cell} python
sd = np.array([0.1, 0.05, 0.25])
```

The procedure is:

```{code-cell} python
D = np.identity(3)
np.fill_diagonal(D, sd)
cov = np.matmul(np.matmul(D, corr), D)
```

The covariance matrix is:

```{code-cell} python
print(cov)
```

Let's simulate a bunch of multivariate normal data from this matrix.
We will sample from a distribution with means of 0.

```{code-cell} python
data = np.random.multivariate_normal([0, 0, 0], cov, size=10000)
```

The simulated means are:

```{code-cell} python
print(data.mean(axis=0))
```

The simulated covariances are:

```{code-cell} python
print(np.cov(data.T))
```

The simulated correlation coefficients are:

```{code-cell} python
print(np.corrcoef(data.T))
```


