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

(demes_event_timings)=

# Timing of events in demographic models

This document provides several examples showing *when* events happen in demographic models.

We will define very simple setups that allow very simple assertions about their outcomes.
While most simulations will involve more complex setups, the advantage of simplicity here
is that we can easily validate the results because we can easily reason about the models.

These examples come from the `fwdpy11` test suite.

The key concept here is that the [demes specification](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html)
defines an epoch as a half open interval on `[present, past)`.
Our reasoning about the outcomes of data that we record during simulation involves breaking our models up into half open intervals and asking questions about them.

Why do we spend time going over this?
Imagine that you want to know the details about the ancestry of a given generation in your simulation.
For example, the ancestry of the generation born during a pulse.
One way to do that would be to preserve both the pulse generation and its parent generation as ancient samples.
(See [here](recorders_vignette).)
To do such sampling you need to know when to sample!

## Some setup

During some of our simulations, we will record the current time point
and the part of the metadata containing individual parents.

We can use a simple data class to record this information:


```{literalinclude} ../../tests/test_demes_event_timings.py
:language: python
:lines: 8-20
```

There will be no genetics happening during the simulation.
All that happens is that we apply our recorder.
We will always "burn in" the demographic model with 100
generations of matings.

```{literalinclude} ../../tests/test_demes_event_timings.py
:language: python
:lines: 24-42
```

```{code-cell} python
---
tags: ['remove-input']
---
import demes
import demesdraw
```

## Intervals when demes exist

In the following graph, demes completely replace one another over time:

```{literalinclude} ../../tests/demes_event_examples/deme_existence.yaml
:language: yaml
```

The following graphic depicts the model:

```{code-cell} python
---
tags: ['remove-input']
---
graph = demes.load("../../tests/demes_event_examples/deme_existence.yaml")
demesdraw.tubes(graph);
```

The intervals of a deme's existence are from $[i, j)$ generations ago, where $j > i$.
Therefore `deme2` exists from $[0, 25)$ generations ago, etc..

During simulation, a node's birth time and deme are recorded.
Therefore, we can export the simulation to a tree sequence and check that the node
times match up with their expected deme labels:

```{literalinclude} ../../tests/test_demes_event_timings.py
:language: python
:lines: 83-93
```

## A single pulse event

The following graph defines a two deme model with a single pulse event 10 generations ago.


```{literalinclude} ../../tests/demes_event_examples/single_pulse.yaml
:language: yaml
```

The following graphic depicts the model:


```{code-cell} python
---
tags: ['remove-input']
---
graph = demes.load("../../tests/demes_event_examples/single_pulse.yaml")
demesdraw.tubes(graph);
```

Given this model:

* The time interval $[0, 10)$ generations ago is the *post* pulse interval 
  (thinking forwards in time).
  During this interval, all individuals in `deme0` have parents from that same deme.
  Likewise for `deme1`.
* All individuals born in the pulse generation (10 generations ago) have all of their
  parents drawn from `deme0`.
  The reason for this is that there are no events affecting the ancestry of `deme0` and
  the proportion of `deme1` ancestry that comes from `deme0` is 1.0, which is all of it.
* The time interval $[11, \infty)$ generations ago is the *pre* pulse interval 
  (thinking forwards in time).
  During this interval, all individuals in `deme0` have parents from that same dime.
  Likewise for `deme1`.

Using the reasoning about the pre and post epochs with respect to the timing of 
the pulse, we can make the following assertions:

```{literalinclude} ../../tests/test_demes_event_timings.py
:language: python
:lines: 45-61
```

Let's work through the logic behind assertions in detail:

* The size of each deme is 100 diploids.
* The metadata are ordered by deme id.
  Thus, the first 100 entries are the indexes parents of individuals in `deme0`, etc..
* Due to the metadata ordering, an index $< 100$ (or, $[0, 100)$) indicates
  a parent from `deme0`.
  An index in $[100, 200)$ indicates a parent from `deme1`.


## A multi-generation "burst" of migration

This model is a slight twist on the previous.
Rather than a single generation pulse, there are several consecutive generations of migration
in one direction:

```{literalinclude} ../../tests/demes_event_examples/burst_of_migration.yaml
:language: yaml
```
The following graphic depicts the model:

```{code-cell} python
---
tags: ['remove-input']
---
graph = demes.load("../../tests/demes_event_examples/burst_of_migration.yaml")
demesdraw.tubes(graph);
```

In this model, migration occurs during the interval $[5, 10)$ generations ago.
Because the migration is one way (`deme0` to `deme1`) and the rate is 1.0, the parents
of all offspring born in this interval have their parents sampled from `deme0`.

Outside of this interval, the parents of any individual are drawn from that
individual's deme.

The following code asserts that these intervals are what we expect:

```{literalinclude} ../../tests/test_demes_event_timings.py
:language: python
:lines: 64-80
```
