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

(demes_vignette)=

# Using `demes` to specify demographic models.

:::{note}
The lists of objects described here are passed to the `demography` and `simlen` parameters when initializing instances of {class}`fwdpy11.ModelParams`.
:::

Starting with version `0.14.0`, it is possible to specify a demographic model using the [`demes`](https://popsim-consortium.github.io/demes-docs/main/) specification.
This method of specifying models uses a `YAML` syntax, which is described in detail in the [`demes specification`](https://popsim-consortium.github.io/demes-spec-docs/main/reference.html).
You may find several examples in the [gallery](https://popsim-consortium.github.io/demes-spec-docs/main/gallery.html#sec-gallery).

There are several advantages to using `demes`:

* The `YAML` files can be used by many pieces of related software, including [`msprime`](https://tskit.dev/msprime/docs/stable/) and [`moments`](https://moments.readthedocs.io).
* The specification is simpler than the low-level `API` provided by any of these tools.
* You'll get an extra layer of model validation when `demes` loads your model prior to it being converted into `fwpdy11` objects.
* Tools to visualize the models are under active development.

## YAML file input

The following `YAML` specifies the human out-of-Africa model from {cite}`Gutenkunst2009-wd`:

```{literalinclude} gutenkunst_ooa.yml
:language: yaml
```

We can generate demographic models directly from these `YAML` files using {func}`fwdpy11.discrete_demography.from_demes`, which handles the conversion to the low-level objects described {ref}`here <softselection>`:

```{code-cell} python
import fwdpy11

model = fwdpy11.discrete_demography.from_demes("gutenkunst_ooa.yml")
```

:::{note}
Be sure to read the documentation for {func}`fwdpy11.discrete_demography.from_demes`!
There are important options concerning the run time of the simulation, etc.
:::

The return value is an instance of {class}`fwdpy11.demographic_models.DemographicModelDetails` and may be passed as the `demography` keyword argument to initialize an instance of {class}`fwdpy11.ModelParams`.
To extract the simulation length to generate the `simlen` parameter of your {class}`fwdpy11.ModelParams` instance:

```{code-cell} python
model.metadata["total_simulation_length"]
```

If we print the `model` object, we will see how much logic we'd have had to use in order to implement this model using the low-level object `API`:

```{code-cell} python
print(model.asblack())
```

It is hopefully clear how much simpler it is to use the `demes` `YAML` specification!

## Working with graphs

You may also build models by creating {class}`demes.Builder` objects using the `demes` `API`.
In general, we feel that the `YAML` method will be less error-prone and therefore the preferred approach.
But, for example:

```{code-cell}
import demes

builder = demes.Builder(description="test demography", time_units="generations")
builder.add_deme(
    name="deme",
    epochs=[
        dict(start_size=1000, end_time=100),
        dict(start_size=2000, end_time=0),
    ],
)

graph = builder.resolve()

model = fwdpy11.discrete_demography.from_demes(graph)

print(model.asblack())
```

Again, it is simpler to build up the demography using `demes` than it is using the `fwdpy11` objects directly.
