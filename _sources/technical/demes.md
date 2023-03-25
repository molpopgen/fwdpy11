(demes_details)=

# Details of `demes` integration

This page gives various technical details regarding `fwdpy11` and [`demes`](https://popsim-consortium.github.io/demes-docs/main/) models.
See {ref}`the vignette <demes_vignette>` for an overview.

## "Burn-in" times

When specifying a burn-in time for a population, that time
corresponds to epochs in the past where the demographic history
of all ancestral demes (those with start times of infinity) are constant.

For example, let's revisit the model from {ref}`before <demes_vignette>`:

```yaml
description: An example demes model
time_units: generations
demes:
 - name: ancestor1
   epochs:
    - start_size: 100
      end_time: 50
 - name: ancestor2
   epochs:
    - start_size: 250
      end_time: 50
 - name: admixed
   start_time: 50
   ancestors: [ancestor1, ancestor2]
   proportions: [0.90, 0.10]
   epochs:
    - start_size: 100
```

Here, we have two ancestral demes, each with a single epoch.
In other words, the size history of these two demes is constant
from time 50 in the past back to "infinity".

For cases like this one, any "burn in" time is simply the history of the
ancestral demes.

For such a model, 49 generations ago will correspond to first time offspring
are born after the burning phase has ended.

Let's modify the model to add a second epoch to ``ancestor2``:

```yaml
description: An example demes model
time_units: generations
demes:
 - name: ancestor1
   epochs:
    - start_size: 100
      end_time: 50
 - name: ancestor2
   epochs:
    - start_size: 100
      end_time: 60
    - start_size: 250
      end_time: 50
 - name: admixed
   start_time: 50
   ancestors: [ancestor1, ancestor2]
   proportions: [0.90, 0.10]
   epochs:
    - start_size: 100
```

Now we have a size change in the ancestral population 60 generations ago.
That event is relevant to contemporary patterns of variation.
Therefore, our burn-in phase must end one generation prior to the size change.

See {func}`fwdpy11.ForwardDemesGraph.from_demes` for details on how to set
the burn-in times.
The {ref}`vignette <demes_vignette>` has examples.

## Selfing

``fwdpy11`` simulates according to Wright-Fisher dynamics.
This means that parents are chosen with replacement and proportionally to their fitnesses.

When an epoch has a selfing rate of zero, sampling with replacement actually implies
a residual selfing rate of {math}`O(1/N)` (assuming that individual fitnesses are 
all close to 1.0).

Further, an epoch with a selfing rate {math}`0 < S \leq 1/N` will have a realized
selfing rate *larger than* {math}`S`.
(An individual will not self with probability {math}`1-S` and will then be chosen
again as a parent with probability {math}`O(1/N)`.)

If this behavior is not desired, pass `allow_residual_selfing = False` when constructing
instances of {class}`fwdpy11.ModelParams`.
Passing ``False`` will prevent an individual from being picked twice as a parent.
See the documentation of {class}`fwdpy11.ModelParams` for details about errors
that may occur if deme sizes get down to a single individual.
