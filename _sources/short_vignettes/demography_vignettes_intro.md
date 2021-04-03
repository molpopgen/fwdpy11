(demography_vignettes_intro)=

# Introduction

## Discrete demes and Wright-Fisher life cycles

`fwdpy11` supports two methods to define demographic models involving discrete demes.
The first method uses Python code and a low-level, object-oriented `API`.
You will find detailed documentation for this approach {ref}`here <softselection>`.

A simpler method is to specify the demographic model via a `YAML` file and build a model using {func}`fwdpy11.DiscreteDemography.from_demes`.
This approach is described in a later vignette and is the recommended method because it will be the least error-prone.

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
