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

(demography_debugger_vignette)=

# Debugging demographic models

You may use {class}`fwdpy11.DemographyDebugger` to check the validity of demographic models.
The following code block applies this class to the model we built using `demes` in {ref}`an earlier vignette <demes_vignette>`. 

For the debugger to do its job, it needs to know the initial sizes of each deme.
For models generated via {func}`fwdpy11.discrete_demography.from_demes`, the meta data contains a {class}`dict` mapping the integer `ID` for each deme to its initial size.
We use that `dict` below to get the required list of initial_sizes.

```{code-cell}
import demes
import fwdpy11

builder = demes.Builder(description="test demography", time_units="generations")
builder.add_deme(
    name="deme",
    epochs=[
        dict(start_size=1000, end_time=100),
        dict(start_size=2000, end_time=0),
    ],
)

graph = builder.resolve()

demog = fwdpy11.discrete_demography.from_demes(graph)

initial_sizes = [demog.metadata["initial_sizes"][i] for i in sorted(demog.metadata["initial_sizes"].keys())]

dbg = fwdpy11.DemographyDebugger(initial_sizes, demog.model)

```

The `dbg` variable also contains a "report" describing the model:

```{code-cell} python
print(dbg.report)
```


