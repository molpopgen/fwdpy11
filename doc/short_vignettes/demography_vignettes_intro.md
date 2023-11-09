(demography_vignettes_intro)=

# Introduction

## Discrete demes and Wright-Fisher life cycles

Demographic models of discrete demes are declared according to the
[demes specification](https://popsim-consortium.github.io/demes-spec-docs/main/tutorial.html).
See {class}`fwdpy11.ForwardDemesGraph` and {func}`fwdpy11.discrete_demography.from_demes` for details.

Please cite the `demes` paper {cite}`Gower2022-xx` when using `fwdpy11`.

## A note of caution

Implementing demographic models is difficult.
There are two sources of difficultly.
The most practical of the two is the software interface.
For example, Hudson's `ms` {cite}`Hudson2002-oo` requires manual specification of all events on the command line.
Complex models required command lines that potentially spanned multiple screens.
(See {cite}`Hernandez2008-cc` for a forward‐time simulation with a qualitatively similar interface.)
Object-oriented interfaces, such as that provided by `msprime` {cite}`Kelleher2016-cb`, can be more readable, yet they will also spin out of control into pages of code for complex models.
These practical issues have had real consequences, leading to errors in the literature {cite}`Ragsdale2020-gl`.

The second difficulty concerns the correctness of the output of a properly-specified model.
By this, we are asking if you are getting the "right" results (in distribution) from your model?
Where possible, we test results from models statistically, using a [separate repository](https://github.com/molpopgen/fwdpy11_statistical_tests).
That repository compares results from demographic simulations to a combination of predictions from analytical theory or numerical calculations performed using [moments](https://moments.readthedocs.io) {cite}`Jouganous2017-tg`.
