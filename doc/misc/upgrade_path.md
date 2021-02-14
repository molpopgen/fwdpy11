(upgrade_path)=

# Upgrade path

This document outlines how to upgrade existing scripts to new versions of fwdpy11.  This guide is likely
imperfect/incomplete.

## 0.12.0

This version removes simulation without tree sequence recording.
The fix is to use tree sequence recording.

## 0.8.2

* The first `kwarg`/positional argument for initializing a
  {class}`fwdpy11.DemographyDebugger` has been renamed `initial_deme_sizes`.
  As the argument is required, and previously only took one possible type,
  we expect this change to not really break anyone's code.

## 0.8.0

* {class}`fwdpy11.DiscreteDemography` can no longer be initialized with a `numpy` array as a positional
  argument. Now, pass it as the value to the `set_deme_sizes` keyword argument.
* Initialization of {class}`fwdpy11.SetMigrationRates` has changed for the case
  of resetting the entire migration matrix. See {ref}`here <migration>`.
* The `shape` `kwarg` to initialize a {class}`fwdpy11.GammaS` has been
  renamed `shape_parameter`.
* The `matrix` `kwarg` to initialize a {class}`fwdpy11.MultivariateGaussianEffects`
  has been renamed `cov_matrix`.
* The `kwarg` `gw2w` for genetic value objects `init` methods has been replaced with `gvalue_to_fitness`.
* The `when` `kwargs` for {class}`fwdpy11.Optimum` and {class}`fwdpy11.PleiotropicOptima` `init` methods is now the **third** positional argument rather than the first.  This change should be caught at run time (the new classes won't auto-convert `float` to `int`, and is fixed by adding parameter names when initializing instances of these types.
* The type of {attr}`fwdpy11.TableCollection.input_left` and {attr}`fwdpy11.TableCollection.output_right` has changed
  to a simple index vector ({pr}`529`).  This change is a BIG memory savings when simulating large regions.

## 0.5.0

The following functions and types previously required a {class}`fwdpy11.MutationVector` argument, but no longer do:

* {class}`fwdpy11.VariantIterator`
* {class}`fwdpy11.DataMatrixIterator`
* {func}`fwdpy11.data_matrix_from_tables`

The extra argument could be eliminated due to the new attributes added to {class}`fwdpy11.MutationRecord`.

## 0.2.0

This release also separates out the data representing a diploid into two classes, {class}`fwdpy11.DiploidGenotype` and {class}`fwdpy11.DiploidMetadata`. 

This release contains major changes to how genetic values are calculated and to how simulations parameters are stored.
These changes are major *simplifications* to the package. 

The changes to how diploid data are stored completely changes how custom genetic values calculations are implemented.

Another major change is that genetic value and noise functions are no longer allowed to be written in Python.  We may
bring that back in a later release.

`fwdpy11.sampling.DataMatrix` has been completely refactored.  See {ref}`datamatrix` for overview of current API.

The function `fwdpy11.sampling.matrix_to_sample` now returns a tuple with two elements, which represent neutral
and selected gentoypes, respectively.  The previous  API made you choose neutral or selected for the return value, which
was a list.

Support for tree sequences will likely have a big impact on how you think about carrying out simulations.  See {ref}`tsoverview`
and {ref}`ts-data-types` for details.

## 0.1.4

### Changes to DataMatrix

The member types `fwdpy11.sampling.DataMatrix.ndim_neutral` and  `fwdpy11.sampling.DataMatrix.ndim_selected` are now read-only attributes.  In previous versions, they were functions.  To upgrade, simply remove any trailing `()`. In other words change this:

```{code-block} python

x.ndim_neutral()

```

To this:

```{code-block} python

x.ndim_neutral

```

The properties `fwdpy11.sampling.DataMatrix.neutral` and `fwdpy11.sampling.DataMatrix.selected` are now
writeable.  This allows you to recode the data as needed.  For example, if you wish to swap the 0/1 values for a column,
subtract 1 then multiply by -1.  The result will affect the data stored on the C++ side.


