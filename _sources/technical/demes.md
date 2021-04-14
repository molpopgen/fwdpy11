(demes_details)=

# Details of `demes` integration

This page gives various technical details regarding `fwdpy11` treats [`demes`](https://popsim-consortium.github.io/demes-docs/main/) models.
See {ref}`the vignette <demes_vignette>` for an overview.

## Implementation of mass migration events (a.k.a. pulse migrations)

The low-level discrete demography `API` provides {class}`fwdpy11.MassMigration` to implement single-generation bulk movements of lineages.
See {ref}`here <massmigrations>` for details.

For models generated with `demes`, we instead use a single-generation change to the migration matrix.
It turns out that this method greatly simplifies the implementation of {func}`fwdpy11.discrete_demography.from_demes`.

While biologically and mathematically equivalent, these two approaches would give different results for the same random number seed.
