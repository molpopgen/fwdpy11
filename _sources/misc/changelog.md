# Changelog

Major changes are listed below.  Each release likely contains fiddling with back-end code,
updates to latest `fwdpp` version, etc.

## 0.21.3

Deprecations

* Constructing {class}`fwdpy11.ModelParams` without a demographic object is now deprecated
  and raises a warning.
  This will be a hard error in a future release.
  PR {pr}`1219`

Documentation

* Minor updates to manual.
  Some formatting issues fixed.
  PR {pr}`1222`
  PR {pr}`1223`

Python back end

* Use `tskit.NODE_IS_SAMPLE` constant instead of relying on the numeric value
  when exporting to `tskit`.
  PR {pr}`1218`

## 0.21.2

Bug fixes

* Fix internal error validating distributions of effect sizes for multi-deme models
  with correlations in effect sizes between demes.
  PR {pr}`1210`.

Examples

* Update examples to latest API and run them in CI.
  PR {pr}`1211`.

Dependencies

* Bump pybind11 to 2.11.11.
  PR {pr}`1172`.
* Require tskit >= 0.5.6
  PR {pr}`1206`.
* Add cmake dependency to pyproject.toml.
  PR {pr}`1177`.
* Bump `demes-forward-capi` to depend on latest version.
  PR {pr}`1182`.

Deprecations

The following are deprecated:

* `fwdpy11.GSS`
* `fwdpy11.GSSmo`
* `fwdpy11.MutivariateGSS`
* `fwdpy11.MutivariateGSSmo`

Their functionality is replaced with {class}`fwdpy11.GaussianStabilizingSelection`.

PR {pr}`1166`.
Issue {issue}{`463`}.

## 0.21.0

Breaking changes

* Refactor {class}`fwdpy11.conditional_models.SimulationStatus` as an enum.
  This change makes correct use much easier.
  PR {pr}`1161`.
  Issue {issue}`934`.

Fixes

* {class}`fwdpy11.conditional_models.ConditionalModelOutput` now contains
  fields allowing one to distinguish if the tracked mutation is still present
  in {attr}`fwdpy11.DiploidPopulation.mutations` and/or {attr}`fwdpy11.DiploidPopulation.fixations`
  PR {pr}`1163`.
  Issue {issue}`1160`.

New features

* {func}`fwdpy11.DiploidPopulation.create_from_tskit` is now able to restore
  individual metadata, populating {attr}`fwdpy11.DiploidPopulation.diploid_metadata` 
  and {attr}`fwdpy11.DiploidPopulation.ancient_sample_metadata`. 
  PR {pr}`1157`.
  Issue {issue}`1130`.

CI changes

* Remove tests using `conda` environments.
  This removal was a pragmatic decision to speed up CI time.
* Test macOS/x86 using brew instead of `conda`
* Reduce number of work flows run for PRs into the `dev` branch.
  Several work flows are only needed when merging into `main`.
* Reduce number of work flows running upon push to `main`/`dev`.
  We now use branch protection, so in theory any changes merging
  have passed their CI requirements.

## 0.20.1

Bug fixes

* Fix attribute errors when validating recombination intervals.
  PR {pr}`1155`
  Issue {issue}`1153`

Build system and CI.

* Update pip usage and ensure reproducible builds of rust code.
  PR {pr}`1156`

## 0.20.0

Documentation

* Update docs on user environments.
  PR {pr}`1148`

Back end

* Disable LTO during CI
  PR {pr}`1147`
* Improve rust handling during CI
  PR {pr}`1146`
* Bump `demes-forward-capi` dependency.
  PR {pr}`1145`
* Bump `demes-spec` submodule.
  PR {pr}`1144`

## 0.20.0a1

New features

* Add {class}`fwdpy11.BinomialIntervalMap`.
  PR {pr}`1142`

Back end

* Improve test suite run times.
  PR {pr}`1136`
* Refactor some C++ test fixtures.
  PR {pr}`1139`
* Python tests now treat warnings as errors.
  PR {pr}`1135`

Documentation

* Improve docs related to {class}`fwdpy11.ForwardDemesGraph`.
  PR {pr}`1133`
  PR {pr}`1137`

## 0.20.0a0

New features

* {func}`fwdpy11.DiploidPopulation.create_from_tskit` now uses the `generation` field from top-level metadata to restore all times to their forward-in-time values.
  PR {pr}`1129`
* Discrete demographic models are now managed by
  {class}`fwdpy11.ForwardDemesGraph`
  See {ref}`here <demes_vignette>`.
  PR {pr}`1101`
  PR {pr}`1103`
  PR {pr}`1121`
  PR {pr}`1123`
  PR {pr}`1124`
* Add support for Python 3.11
  PR {pr}`1108`
* Allow suppressing of residual selfing.
  See {ref}`here <demes_details>`.
  PR {pr}`1117`
  PR {pr}`1118`

Deprecations

* Deprecate use of {class}`fwdpy11.Region` for recombination
  PR {pr}`1126`

Breaking changes

* Drop support for Python 3.7
  PR {pr}`1108`
* Removed deprecated discrete demography API
  PR {pr}`1101`
  PR {pr}`1116`

Documentation

* Give link to examples of lists of "gvalue" objects
  PR {pr}`1127`
* Rewrite `demes`-related documentation
  PR {pr}`1124` 

## 0.19.10

Documentation

* Rewrite {ref}`import_mutations_from_tskit_vignette`.
  PR {pr}`1119`
  Issue {issue}`1104`

## 0.19.9

Big fix

* Use the mutation's time field rather than
  the metadata "origin" time when creating
  a DiploidPopulation from a tskit tree sequence.
  PR {pr}`1113`

## 0.19.8

Big fix/new feature (take your pick)

* {func}`fwdpy11.DiploidPopulation.create_from_tskit` now
  handles tree sequences were all mutations come from
  adding them "manually" in tskit
  (see {ref}`import_mutations_from_tskit_vignette`)
  **OR** from a previous run of fwdpy11.
  PR {pr}`1110`
  Issue {issue}`1109`

## 0.19.7

Bug fix

* Fix bug in ModelParams processing of non-Poisson
  crossover regions.
  PR {pr}`1106`
  Issue {issue}`1105`

## 0.19.6

Bug fix

* Fix error in conditional models when sampling ancient
  samples throughout the focal mutation sojourn.
  PR {pr}`1094`
  Issue {issue}`1093`

## 0.19.5

Bug fix

* Fix error in BinomialInterval and BinomialPoint where
  they always returned a breakpoint.
  This bug was introduced in 0.19.4.
  PR {pr}`1088`

## 0.19.4

Performance

* Greatly improve the performance of recombination maps.
  PR {pr}`1087`.
  Issue {issue}`1082`.

Python API

* Add :class:`fwdpy11.ForwardDemesGraph`.
  PR {pr}`1077`.

Back end changes

* Develop infrastructure to use demes models directly for evolving populations.
  PR {pr}`1069`.
  PR {pr}`1070`.
  PR {pr}`1071`.
  PR {pr}`1072`.
  PR {pr}`1073`.
  PR {pr}`1074`.
  PR {pr}`1075`.
  PR {pr}`1076`.

## 0.19.3

Deprecations

* The low-level Python API for defining demographic events
  is now deprecated and will raise warnings.
  No action has to be taken yet.
  There is no replacement API to migrate to.
  PR {pr}`1064`.
    
Back end changes

* Prepare to refactor C++ code to evolve directly from `demes` graphs.
  PR {pr}`1063`.
  PR {pr}`1067`.
* Boolean flags affecting simulation behavior collected into `struct`.
  PR {pr}`1061`.

## 0.19.2

Documentation

* Add link between vignettes regarding
  importing tree sequences from tskit.
  PR {pr}`1054`.
* Force notebooks to be executed during doc builds.
  PR {pr}`1058`.

Testing/CI

* Docs workflow updated to run with Python 3.10.
  PR {pr}`1056`.

## 0.19.1

Dependencies

* Bump jupyter-book to 0.13.1
  PR {pr}`1052`

Deployment

* Fix source dist generation during wheel building
  PR {pr}`1051`

## 0.19.0

Differences from the last stable release also include entries
for all alpha releases listed below

Documentation

* Add new vignettes.
  See {ref}`here <varying_fitness_conditions_across_demes_intro>`
  PR {pr}`950`

Bug fix

* Add runtime checks for invalid models where effect sizes
  vary across demes.
  PR {pr}`1047`

Testing

* Disable wheel build workflow for PRs
  PR {pr}`1048`


## 0.19.0a5

Deployment

* Add missing , to docker tag names
  PR {pr}`1042`

Testing

* Remove stray return statement in a test.
  PR {pr}`1043`


## 0.19.0a4

Deployment

* Add a version name tag to releases
  PR {pr}`1041`

Documentation

* Describe podman and give singularity example in deployment section.
  PR {pr}`1040`

Testing

* Fix logic error in testing mutation import from tskit.
  PR {pr}`1039`

## 0.19.0a3

Deployment

* Update docker base image to `ubuntu:jammy`.
  PR {pr}`1031`
* Docker images now installs into a venv rather than globablly.
  PR {pr}`1031`
* Fix docker tags so that `latest` means the latest release.
  PR {pr}`1036`

Documentations

* Document how to build core library in isolation.
  PR {pr}`1038`
* Document `cbindgen` requirement for development.
  PR {pr}`1030`

Back end changes

* C++ back ends of Region/Sregion now validate their coordinates at run time.
  PR {pr}`1032`
* Add wrapper to `gsl_ran_flat` to avoid returning max value.
  PR {pr}`1035`

Testing

* Separate Python tests requiring compiled C++ modules to 
  another directory.
  PR {pr}`1037`

## 0.19.0a2

Documentation

* Rewrite {ref}`vignette <import_mutations_from_tskit_vignette>`.
  PR {pr}`1026`.

## 0.19.0a1

New features

* Ability to import tree sequences with mutations from tskit.
  A new {ref}`vignette <import_mutations_from_tskit_vignette>` describes the procedure.
  PR {pr}`1026`.

Bug fixes

* Fix bug in sweep models.
  The fix in PR {pr}`1020` was incomplete.
  Tree span is now accounted for.
  PR {pr}`1026`.
* Fix Python-side constructors in Mutation re: integer types.
  No previous results were affected.
  PR {pr}`1023`.

Documentation

* Update developer docs.
  PR {pr}`1021`.

## 0.19.0a0

Bug fixes

* Fix bug in sweep models.
  The previous implementation chose acceptable branches for mutation placement uniformly.
  The new method chooses proportionally to branch length.
  PR {pr}`1020`.

API changes

* Remove deprecated `WrappedTreeSeqeuence`
  PR {pr}`978`.
* Add Python-level checks that "regions" don't contain positions extending past the genome ends.
  PR {pr}`987`.

Dependencies

* Bump pybind11 to 2.10.0
  PR {pr}`986`.

Build system and CI

* Remove all use of automake/autoconf.
  PR {pr}`971`.
* Build core library as shared library.
  PR {pr}`961`.
* Allow for standalone cmake builds.
  PR {pr}`966`.
* Add rust crate `demes-forward-capi` to core lib.
  PR {pr}`962`.
  PR {pr}`1001`.
* Build C++ tests with cmake.
  PR {pr}`967`.
* Use `python -m build .` for Ubuntu CI.
  PR {pr}`970`.
* Build all C++ targets in one command on Ubuntu CI.
  PR {pr}`972`.
* Test docker work flow on CI.
  PR {pr}`979`.
* Add cancel actions to all work flows.
  PR {pr}`982`.
  PR {pr}`985`.
* Wheels are now based on `manylinux_2_28`
  PR {pr}`994`.
* move `fwdpy11/src` to `cpp/`
  PR {pr}`996`.
* Add dependabot for GitHub actions
  PR {pr}`1002`.
* Fix warnings from `python -m build .`
  PR {pr}`1018`.

Back end changes

* Move evolve code to core library.
  PR {pr}`974`.
* Clean up implementation of individual and mutation metadata decoding.
  PR {pr}`1015`.
* Tidy up a lot of `mypy` errors.
  PR {pr}`1016`.

Testing

* Remove many C++ files from tests/
  PR {pr}`976`.
* Test that we can load YAML files from demes-spec
  into our rust-based ForwardDemesGraph.
  PR {pr}`975`.
* move `fwdpy11.ezparams` module to tests/
  PR {pr}`977`.
* Add C++ tests of our new "forward demes graph".
  PR {pr}`980`.
  PR {pr}`983`.
  PR {pr}`984`.
  PR {pr}`989`.
  PR {pr}`999`.

## 0.18.3

Bug fixes:

* Fix error in mapping population names to integer IDs when importing `demes` models.
  Issue {issue}`990`.
  PR {pr}`991`.

## 0.18.2

Back end changes

* Simplify the back-end behind {class}`fwdpy11.MutationDominance`.
  PR {pr}`947`.
* Refactor a class decorator used to pickle/unpickle `attrs`-based classes.
  PR {pr}`948`.
* Reject any `demes` models that involve multiple pulse events into the
  same source deme at the same time.  Such models are better written
  in a different way.
  PR {pr}`953`

Bug fixes

* Fixed edge case when importing demes models with burn-in time of 0.
  PR {pr}`949`.

Dependencies

* Pinned pybind11 to 2.9.1.
  PR {pr}`955`.
* Use [build](https://pypi.org/project/build/) for wheels.
  PR {pr}`958`.
* Bump `tskit` to `=~ 0.5.0`.
  PR {pr}`959`.

Build system and CI

* Fix dependency caching for Ubuntu tests.
  PR {pr}`957`.
* Fix wheel builds for macosx.
  PR {pr}`941`.

## 0.18.1

Bug fixes

* Fix display of `weight` for distribution-of-effect-size classes.
  Issue {issue}`885`. PR {pr}`943`.

## 0.18.0

Breaking changes

* {func}`fwdpy11.TableCollection.fs` no longer accepts more than two sample sets.
  This change allowed us to drop `sparse` as a dependency that was causing
  headaches when new Python point releases come out.
  PR {pr}`924`. Issues {issue}`876`, {issue}`919`.
* {func}`fwdpy11.evolvets`: the default for `suppress_table_indexing`
  has changed from `False` to `True`.
  This change will generally result in faster simulations, but will
  break work flows that relied on accessing the trees during simulation.
* Removed deprecated attributes from distribution of effect size types.
  These attributes have been undocumented for so long that no one's code
  should break.
  PR {pr}`938`. Issue {issue}`886`.

Behavior changes

* `fwdpy11.DemographyDebugger.report` raises a warning
  to state that it has not been implemented.
  It returns a string with a message to that effect.

Bug fixes

* Fixed a bug in event time handling from `demes` models.
  This bug did not lead to incorrect results as the C++ back end caught the
  invalid models at run time.
  PR {pr}`891`.
  Issue {issue}`881`.
  {user}`apragsdale`,
  {user}`molpopgen`.
* Fix lambda capture in `examples/plugin/gvalue_recorder.cc`.
  The previous code was not compatible with current versions of `pybind11`.
  PR {pr}`921`. Issue {issue}`920`.
* Fixed handling of non-integer times in `demes` models.
  PR {pr}`930`. Issue {issue}`929`.

User interface improvements

* Region coordinates are now validated before starting the simulation.
  PR {pr}`909`.
  Issue {issue}`908`.

Testing

* Added many more tests related to models defined using `demes`
  PR {pr}`891`.
  PR {pr}`931`.
  Issue {issue}`890`.
  {user}`apragsdale`,
  {user}`molpopgen`.

Back end changes

* Add new Python class intended to help validate demographic models and
  aid in importing models from `demes`.
  PR {pr}`900`.
  PR {pr}`904`.
* Completely rebuild {class}`fwdpy11.DemographyDebugger`.
  PR {pr}`906`.
* Added infrastructure to separate out back end from Python-specific C++ code.
  PR {pr}`936`.
* Assert that epochs from `demes` models are at least 1 generation long.
  PR {pr}`931`
* Added infrastructure to separate out back end from Python-specific C++ code.
  PR {pr}`936`.

Dependencies

* Removed use of `pip-tools` in favor of a single `requirements/development.txt` file with loose pinning.
  Issue {issue}`877`, PR {pr}`896`
* Bump minimum `pybind11` to 2.9.0.
  PR {pr}`922`.

Build System

* `setup.cfg` now pins minimum and maximum Python versions.
  PR {pr}`923`.  Issue {issue}`914`.

Deployment

* Update Docker work flow for building wheels to correctly locate requirements file 
  and build for Python 3.10
  PR {pr}`927`
* Update GitHub actions to build wheels for Python 3.10.
  PR {pr}`928`

## 0.17.1

Bug fixes

* Fix bug in fixation checking in `fwdpy11.conditional_models`.
  Previous results were not incorrect.
  Rather, they took too long to be obtained because some simulations did not terminate early enough.
  Issue {issue}`893`
  PR {pr}`894`

  Miscellaneous

- Improved exception message for invalid migrations during a simulation.
  PR {pr}`892`

## 0.17.0

New release!
This release contains everything from the alpha release series plus the changes noted below.

Documentation:

* Add basic documentation for `fwdpy11.conditional_models`.
  PR {pr}`878`.

Breaking changes with respect to previous alpha releases:

* Change attribute names in object returned by conditional models.
  PR {pr}`878`.

Dependencies

* pin `numpy` in `setup.cfg`.
  This change fixes issues with docker images.
  PR {pr}`880` Issue {issue}`879`

## 0.17.0a4

Dependencies:

* Clean up pinning of `demes` and `demesdraw`.
  PR {pr}`869` 
* Update to `demes` ~=0.2.0.
  PR {pr}`873`
* Update to `tskit` ~=0.4.0
  PR {pr}`874`

## 0.17.0a3

Changes to `fwdpy11.conditional_models`:

* Return evolved instance of {class}`fwdpy11.ModelParams`.
  PR {pr}`866`

Bug fixes

* Fixed bug in updating `fwdpy11.MultivariateGSSmo`.
  PR {pr}`867`.

Back end changes:

* `fwdpy11.GSSmo` and `fwdpy11.MultivariateGSSmo` now handle cases where population start time is greater than zero.
  PR {pr}`867`.

Build system:

* Update requirements.txt and doc/requirements.txt.
  PR {pr}`864`.

Documentation:

* Manual is now built with "nitpick" enabled.
  PR {pr}`865`.
* Updated {ref}`precapitation` to latest `msprime` API.
  Moved "recapitation" concept into a {ref}`shorter vignette <recapitation>`.
  PR {pr}`868`.

## 0.17.0a2

Changes to `fwdpy11.conditional_models`:

* Change name of `kwarg` to `track_added_mutation`.
  PR {pr}`858`

Back end changes:

* Remove extra copies of tskit table collections from the implementation of {func}`fwdpy11.tskit_tools.iterate_timepoints_with_individuals`.
  PR {pr}`859`
  Issue {issue}`851`

Bug fixes:

* Fixed a "use after move" error in C++ code used to pickle `fwdpy11.DiscreteDemograpy`.
  It looks like this bug was introduced back in {pr}`791`, which was part of the 0.16.0 release.
  PR {pr}`857`


## 0.17.0a1

Changes to `fwdpy11.conditional_models`:

* Rename `track_mutation` to `track_added_mutation`.
  PR {pr}`856`

## 0.17.0a0

Bug fixes

* Fix error initializing the founder genome of {class}`fwdpy11.DiploidPopulation`
  Issue {issue}`836`
  PR {pr}`838`
* Fix bug in handling mutation counts when final generation was recorded as "ancient samples".
  This bug resulted in an exception being raised.
  Thus, previous exception-free results were not affected.
  Issue {issue}`844`
  PR {pr}`845`
  PR {pr}`847`

New features

* Add module {mod}`fwdpy11.conditional_models`.
  PR {pr}`828`

Deprecations

* Deprecate `fwdpy11.tskit_tools.WrappedTreeSequence`
  PR {pr}`841`
  PR {pr}`848`

Back end changes:

* Add more runtime checks to {func}`fwdpy11.DiploidPopulation.add_mutation`
  PR {pr}`837`
* Remove unnecessary copying of table collections in {mod}`fwdpy11.tskit_tools`.
  PR {pr}`842`
* Move back-end for ancient sample recording to the population classes.
  Add `fwdpy11.DiploidPopulation._record_ancient_samples`.
  {pr}`853`

## 0.16.2

Documentation

* Document that virtual envs should upgrade `pip`.
  PR {pr}`834`
  Issue {issue}`833`

Packaging

* Fix import of functions in `fwdpy11.demographic_models.human`.
  PR {pr}`835`
  Issue {issue}`832`

## 0.16.1

New features:

* `Mutation.__str__` is now more informative.
  PR {pr}`825`

## 0.16.0

Bug fixes

* An integer type in the infinitely-many sites mutation model was changed from unsigned to signed.
  This does not affect previous results because unsigned overflow doing the "right thing" ended up with final values being correct.
  PR {pr}`766`
  Issue {issue}`765`
* Fix a bug where stopping/restarting the evolution of demographic models at time points
  where a deme goes extinct.
  It is not possible that this bug affected results from earlier versions, as attempting to stop/start at these time points raised exceptions.
  Issue {issue}`775`
  PR {pr}`774`
* Fix bugs in C++ back-end for discrete demographic models.
  In some cases, we were using the wrong vector of deme sizes to update the model,
  leading to runtime exceptions.
  PR {pr}`802`
  PR {pr}`803`
* Fix error in `demes` models where "replacement" models had 1 generation of overlap between ancestral/derived demes.
  Issue {issue}`814`
  PR {pr}`815`
  {user}`apragsdale`
  {user}`molpopgen`

Behavior changes

* If a demographic model is evolved, pickled, unpickled, and then used to evolve,
  it is now possible that exceptions will raise.
  This change is due to the fix for Issue {issue}`775` introduced in PR {pr}`774`.
  See issue {issue}`777` for more background.
* Mass migration events implemented via `fwdpy11.copy_individuals`
  and `fwdpy11.move_individuals` now occur *after* sampling within a generation.  
  This change makes the timings consistent with all other events and also makes
  certain operations easier/feasible.
  {pr}`809`
* Calling {func}`fwdpy11.infinite_sites` during a simulation now raises `RuntimeError`.
  {pr}`820`
  {issue}`769`
* Models imported from `demes` now start the forward-time portion of the model 1 (one) generation before the most ancient end time of an ancestral deme.
  {pr}`818`
  {user}`apragsdale`
  {user}`molpopgen`

New features

* Add {func}`fwdpy11.DiploidPopulation.add_mutation`.
  PR {pr}`764`
  PR {pr}`799`
* Add {class}`fwdpy11.NewMutationData`.
  PR {pr}`764`
* Add `__copy__` and `__deepcopy__` to {class}`fwdpy11.DiploidPopulation`.
  PR {pr}`770`
* Add `__deepcopy__` to `fwdpy11.DiscreteDemograpy`.
  PR {pr}`773`

C++ back-end

* A population can now be checked that it is- or is not- being simulated.
  PR {pr}`762`
* `fwdpy11.discrete_demography.DiscreteDemography` now stores the migration matrix as a stack-allocated object and not a `unique_ptr`.
  PR {pr}`785`
  {issue}`781`

Build system

* All GCC builds and CI tests on Ubuntu + GCC now apply a much stricter set of compiler options.
  {pr}`779` {issue}`778`

Dependencies

* Bump `pillow` version in doc/requirements.txt.
  {pr}`763`
  {pr}`811`
* Bump all lib dependencies and some doc dependencies.
  {pr}`807`

## 0.15.2

Point release

Bug fixes

* `fwdpy11.SetMigrationRates` now uses a tolerance when checking that rates sum to 0 or 1.
  {issue}`787`
  {user}`apragsdale`
  {user}`molpopgen`
* Fix bug where the number of rows in a {class}`tskit.PopulationTable` were incorrect upon export via {meth}`fwdpy11.DiploidPopulation.dump_tables_to_tskit`.
  Issue {issue}`792`.
  Fixed in PR {pr}`793`.

Minor changes

* Remove use of deprecated `numpy` `dtype` in tests.
  Issue {issue}`789`.
  Fixed in PR {pr}`794`
* Added more tests of `demes`-generated models with symmetric migration and individual demes going extinct.
  Issue {issue}`758`
  Closed by {pr}`797`

Dependencies

* Pinned `demesdraw` in `doc/requirements.txt`
  Issue {issue}`790`.
  Fixed in PR {pr}`795`.
* Pin `demes` to `== 0.1.2` in `setup.cfg`.
  Issue {issue}`796`.
  Fixed in PR {pr}`795`.

## 0.15.1

Point release

Bug fixes

* Fix error in decoding provenance rows when initializing `fwdpy11.tskit_tools.WrappedTreeSequence`.
  {pr}`760`

Dependencies:

* `demes` bumped to `>=0.1.2` in requirements files.
  This change is for `demesdraw` and the manual.
  It is still `>= 0.1.1` in setup.cfg.
  {pr}`761`

## 0.15.0

This release is mostly about tskit.
All changes and issues are collected under the `0.15.0` milestone on `GitHub`.

Breaking changes:

* Dropped support for Python 3.6.
  Now, we support 3.7, 3.8, and 3.9.
  {pr}`735`

New documentation:

* Rewrote {ref}`vignette <tskitconvert_vignette>` on exporting data to `tskit`.
* Added {ref}`new vignette <tskit_metadata_vignette>` on working with data exported to `tskit`

These docs were added in {pr}`745`.

New features:

* Added `fwdpy11.tskit_tools.WrappedTreeSequence`.
  {pr}`743`
  {pr}`747`
* {func}`fwdpy11.DiploidPopulation.dump_tables_to_tskit` may now return a `WrappedTreeSequence`.
  {pr}`748`
* Allow row/slice semantics for decoding `tskit` metadata.
  {pr}`734`
* Top-level metadata for `tskit` objects is now filled.
  See {ref}`here <tskitconvert_vignette>` for details.
* Custom deme names can now be added to `tskit` population tables.
  {pr}`742`

Fixes:

* Tables are now indexed after calling {func}`fwdpy11.DiploidPopulation.load_from_file`.
  {pr}`739`
* Various fields of {class}`fwdpy11.tskit_tools.DiploidMetadata` are now populated as `bool` as documented.
  {pr}`742`
* Mutation metadata for `tskit` changed so that exported tree sequences are compatible with `msprime.sim_mutations`.
  {pr}`731`.

## 0.14.1

This is a point release adding more documentation:

* {ref}`Demes vignette <demes_vignette>` updated.
* {func}`fwdpy11.TableCollection.fs` docstring updated regarding some perhaps unexpected behavior of `sparse.COO`.

## 0.14.0

In addition to the changes listed below, several documentation, CI, and deployment changes also happened.
These are collected under the `0.14.0` milestone on `GitHub`.

New features:

* Initial support for demographic models using the `demes` spec.
  See {ref}`here <demes_vignette>`.
  {pr}`710`
  {pr}`712`
  {pr}`713`
  {user}`apragsdale`
* {func}`fwdpy11.DiploidPopulation.dump_tables_to_tskit` has a new option to "destructively" dump tables, which may save a lot of memory.
  {pr}`695`

Python changes:

* Improved support for type hints
  {pr}`690`
  {pr}`692`
  {pr}`694`
  {pr}`696`
* Deprecated `popsizes` and `pself` keyword arguments removed from {class}`fwdpy11.ModelParams`.
  {pr}`703`
* Fixed a design issue in how {class}`fwdpy11.DemographyDebugger` parsed input event lists.
  {pr}`693`

C++ back-end changes:

* Several C++ source files were changed to no longer include `pybind11` headers.
  {pr}`705`
  {pr}`704`
  {pr}`702`
* Removed unused C++ header files
  {pr}`689`

## 0.13.2

Point release

* Small fix to wheel building action.
  {commit}`c31420aede4180becfe2a28936469ed471c6ea41`

## 0.13.1

Point release

* Attempt to fix deployment of various deliverables upon release.
  {pr}`682`

## 0.13.0

API changes

* The `demesizes` keyword is no longer accepted to initialize instances of {class}`fwdpy11.DiploidPopulation`. 
  {pr}`676`
* {meth}`fwdpy11.TreeIterator.nodes` no longer takes an argument.
  {pr}`678`
* The type of {attr}`fwdpy11.Mutation.g` changed from an unsigned integer to a signed integer.
  This change has no practical consequence to user code written in Python.
  {issue}`656`
  {pr}`667`
  {pr}`670`

New features

* The dominance of a mutation may now be specified by a function.
  See {ref}`mutationdominance_vignette`.
  {pr}`590`
  {pr}`629`
  {pr}`630`
* Add {meth}`fwdpy11.TableCollection.build_indexes`.
  {pr}`651`

Improved IDE and editor integration

* Several key classes and functions were reimplemented.
  They are now Python layers on top of the C++.
  These changes allow better type hinting and make the docstrings more easily discoverable, resulting in `jedi`, etc., being better able to use them for auto completion and other tasks.
  {pr}`676`
  {pr}`677`
  {pr}`678`

Packaging changes

* The build system is now PEP 517/518 compliant.
  {commit}`cfd95f57e97d0fff6d4d87d06bd280d2f58ea545`
* `setup.py` no longer imports `pybind11`.
  {pr}`634`

Bug fixes

* Fixed bugs in experimental features involving simulating neutral mutations during simulations recording ancient samples.
  These fixes make these features less experimental and part of the API now.
  {pr}`643`
  {pr}`647`
  {pr}`650`
* Improved input parameter validation for {func}`fwdpy11.evolvets`.
  {pr}`649`
* Fixed a bug where mutations on branches generated by `msprime` had invalid origin times, leading to exceptions when calling {meth}`fwdpy11.DiploidPopulation.dump_tables_to_tskit`.
  {pr}`656`
  {pr}`670`

Deployment changes

* Pushes to `main` and new releases now trigger updates to a Docker image.
  See {ref}`sec_deployment`.
  {pr}`652`
  {pr}`663`
  {pr}`671`
  {pr}`672`
* Binary wheels are automatically deployed for macOS and Linux for new releases.
  {pr}`669`

Documentation changes

* The manual is now built with [jupyter book](https://www.jupyterbook.org). {pr}`653`.
* The manual is now longer hosted on Read The Docs.  
  It has been moved [here](https://molpopgen.github.io/fwdpy11) and is automatically deployed by GitHub actions.
  {pr}`654`.

## 0.12.0

API breakage:

* Simulation without tree sequence recording is no longer supported. {pr}`607`

Bug fixes:

This release fixes two bugs.
Neither bug would affect results.
Rather, they affect the behavior of some advanced features.

* Explicit initialization of a data structure used for simplification is now always done. {issue}`626` {pr}`627`
* Simplification is now allowed to occur after a single generation. {issue}`624` {pr}`627`

New features:

* Recombination breakpoints can be restricted to integer positions. {pr}`612`
* Improved methods for decoding `tskit` meta data. {pr}`609`.

<!---
note:

The meta data decoding is still a moving target/work in progress.
We anticipate further changes as we get a better idea of what works well/is most useful.
-->

This release also includes several updates to dependencies, documentation, and CI testing.
These changes are collected under the `0.12.0` milestone on `GitHub`.

## 0.11.0

Point release.

Behavior changes:

* The events passed to `fwdpy11.DiscreteDemograpy` are now sorted stably
  by time. {pr}`598`

Documentation changes:

* Fixed documentation for `fwdpy11.SetMigrationRates`. {pr}`603`

Back end changes:

* Cleanup some parts of the code base.  {pr}`604`
* Make the C++ back-end of `fwdpy11.MassMigration` have the same initialization
  semantics as other demographic events. {pr}`605`.

## 0.10.1

Point release.

Bug fixes:

* Fix some errors in {class}`fwdpy11.DemographyDebugger`. {pr}`594`

Testing changes:

* We now use GitHub actions instead of Travis. {pr}`593`

Documentation changes:

* We now use `jupyter-sphinx` to execute code in the manual. {pr}`597`

## 0.10.0

Interface changes:

* When dumping tables to `tskit` via {func}`fwdpy11.DiploidPopulation.dump_tables_to_tskit`, we now use the latest metadata methods.
  See {ref}`here <tskit_metadata_vignette>`.
  {pr}`588`
* Most arguments to {func}`fwdpy11.evolvets` are now keyword-only.
  {pr}`584`
* Added {class}`fwdpy11.DiscreteDESD`.
  {pr}`587`

Dependency updates:

* Minimum `tskit` version is now 0.3.2.
* Minimum `attrs` version is now 0.19.2.
* `Sphinx` version pinned to 3.1.
* The warnings for deprecated features got promoted to {class}`FutureWarning`.
  {pr}`583`

C++ back-end changes:

* New method to handle GSL errors. {pr}`574`
* Table collections are now managed by shared pointers. {pr}`582`

## 0.9.0

This release enables custom genetic value models to be implemented in Python.
To do so, the back-end for C++ genetic values was changed in a way that (hopefully!)
future-proofs the API against future changes.  The approach taken to allowing
Python genetic value types evolved quite a bit during development, so we won't
refer to individual pull requests here.  Anyone interested can look at the 0.9.0
milestone on GitHub.

See {ref}`here <gvalues_python>` for the documentation on Python genetic values.

## 0.8.3

* {func}`fwdpy11.DiploidPopulation.dump_tables_to_tskit` now populates
  the provenance table. {pr}`542`
* Improve checking migration rates in {class}`fwdpy11.DemographyDebugger`. {pr}`545`
* {class}`fwdpy11.DemographyDebugger` now makes a deep copy of input. {pr}`546`
* The C++ back-ends for Gaussian stabilizing selection classes got streamlined
  without changing the user interface. {pr}`547`
* Manual got overhauled. {pr}`543`
* Snowdrift example reimplemented using `attrs`. {pr}`548`

## 0.8.2

* Fix issue where {class}`fwdpy11.DemographyDebugger` failed to
  catch populations with empty migration matrix rows after
  mass migration. {pr}`539`
* {class}`fwdpy11.DemographyDebugger` is now implemented
  with `attrs`. {pr}`540`.  This change changes a keyword
  argument for this class.  See {ref}`upgrade guide <upgrade_path>`.

## 0.8.1

* Fixed a back-end bug that could have led to corrupt sample lists for simplification. {pr}`536`.
* Made improvements to memory handling of data structures when simulations end. {pr}`537`.
* Added the three-deme model of Jouganous et al. (2017).
  See `fwdpy11.demographic_models.human.jouganous_three_deme`.
  {pr}`534`

## 0.8.0

Stable release. In addition to what is in the previous alpha releases:

* Memory use is substantially reduced due to some improvements
  in `fwdpp`.  {pr}`533` brings in two changes from `fwdpp`,
  {pr}`molpopgen/fwdpp#287` and {pr}`molpopgen/fwdpp#288`.

This release includes some minor `API` changes.
See the {ref}`upgrade guide <upgrade_path>` for details.

## 0.8.0a1

Second alpha release of 0.8.0:

* Update the `fwdpp` submodule. {pr}`529`
* Update documentation related to genetic maps. {pr}`530`

## 0.8.0a0

This is the first alpha release of 0.8.0.

In addition to what is below, this release contains
a smattering of build system changes, documentation changes,
etc., that are collected under the 0.8.0 milestone on `Github`.

`API`/`UI` changes:

This release brings Python classes that have been reimplemented using [attrs](https://www.attrs.org).  These changes have a lot of benefits:

* A lot of C++ code got removed (yay!) because we can use `attrs` for the pickling
  machinery, `__repr__`, etc..
* We now get *much* nicer `__repr__` for all of the types that get sent into
  instances of {class}`fwdpy11.ModelParams`.

However, these changes required some simplification to the `__init__` methods,
which meant some `API` breakage. See the {ref}`upgrade guide <upgrade_path>`
for details.

This release also removes features deprecated in previous releases. {pr}`482`

Performance improvements:

* Sorting edge tables prior to tree sequence simplification has been replaced
  by an efficient buffering algorithm. {pr}`526`.

New demographic models:

* The [^cite_tennessen2012] model is added via `fwdpy11.demographic_models.human.tennessen`.
  {pr}`479`

Improved behavior:

* Improved warnings about demographic events scheduled to happen
  before the population's current generation. {pr}`495`
* Built-in demographic models now return instances of
  `fwdpy11.demographic_models.DemographicModelDetails`.
  Such instances can be passed as the `demography` keyword argument
  to initialize {class}`fwdpy11.ModelParams`.
  {pr}`509`.
* The "individual" column of a node table is now populated
  when exporting to a {class}`tskit.TableCollection`. {pr}`488`

Changes to implementation of Python classes

* {class}`fwdpy11.ModelParams` has been reimplemented
  using [attrs](https://www.attrs.org). {pr}`484`, {pr}`486`, {pr}`487`.
* Demographic model types are now implemented using [attrs](https://www.attrs.org) and
  inherit from the C++ back-end class. {pr}`492`
* Region types are now implemented using [attrs](https://www.attrs.org) and
  inherit from the C++ back-end class. {pr}`497`
* Genetic value types are now implemented using [attrs](https://www.attrs.org>) and
  inherit from the C++ back-end class. {pr}`504`
* Genetic map unit types are now implemented using [attrs](https://www.attrs.org) and
  inherit from the C++ back-end class. {pr}`506`

C++ back end changes:

* The default C++ language standard is now C++14. {pr}`517`.
* Custom exceptions now have default symbol visibility. {pr}`519`.
* The back-end code for discrete demography got cleaned up. {pr}`521`.
* The `fwdpp` submodule was updated a few times.
  {pr}`489` {pr}`523` {pr}`525`

## 0.7.1

Maintenance release and one new feature:

* Allow the first generation of a simulation to be preserved. PR {pr}`470`
  See {ref}`recapitation`.
* Parameterizing classes like `fwdpy11.GSSmo` is now more Pythonic,
  and some existing `init` methods are deprecated in favor of the
  new approach. PR {pr}`461`.

This release include several other improvements to documentation and user interface.
All changes are backwards-compatible, and deprecation warnings are issued when
necessary.  See the 0.7.1 milestone on `GitHub` for details.

## 0.7.0

Major feature release allowing mutations to have different
effect sizes in different demes.

Bugs fixed:

* Temporal samplers now get the correct offspring metadata in simulations
  with tree sequence recording. {issue}`444`

New features:

* Added {class}`fwdpy11.mvDES`, which allows multivariate distributions of effect sizes
  such that mutations have different effect sizes in different demes. See {ref}`mvdes`
  for details. PR {pr}`443` PR {pr}`452`
* {class}`fwdpy11.GeneticValueToFitnessMap` now records whether or not genetic
  values are mapped to fitness or are a trait value via {attr}`fwdpy11.GeneticValueToFitnessMap.maps_to_fitness`
  and {attr}`fwdpy11.GeneticValueToFitnessMap.maps_to_trait_value`.
  PR {pr}`447`

Other changes (see the 0.7.0 milestone on GitHub)

* This release deprecates several features that are no longer sensible given that most
  simulations will use tree sequence recording.  You will see warnings pop up if you
  use these features (or run the unit tests).  These features will be removed
  in 0.8.0.
* Many back-end changes to the C++ code simplify things in various places.

## 0.6.4

Fixes a bug where the timing of updates to stateful genetic values
was off by one generation:

* {issue}`437`

## 0.6.3

Maintenance release.

This release fixes three bugs. The first two are related to internal
details of book-keeping various data structures:

* {issue}`420`
* {issue}`422`
* {issue}`432`

Other changes:

* `sparse` is added to `install_requires` in `setup.py`.  {issue}`421`
* {class}`fwdpy11.TableCollection`'s validation of genome lengths is improved. PR {pr}`428`
* The C++ base class for a population is now a concrete class rather than a template alias.  This change enables forward declarations in header files. PR {pr}`427`

## 0.6.2

This release changes the migration code to model juvenile migration.
These changes simplify the back end and give the same results (in
distribution).  The relevant PRs are:

* PR {pr}`416`
* PR {pr}`417`

## 0.6.1

This is a maintenance release that clears up a few issues:

* {issue}`246`
* {issue}`280`
* {issue}`339`
* {issue}`365`
* {issue}`386`
* {issue}`397`

The following features are added:

* {attr}`fwdpy11.DataMatrix.neutral_matrix`
* {attr}`fwdpy11.DataMatrix.selected_matrix`
* {func}`fwdpy11.DataMatrix.merge`

## 0.6.0

This is a major feature release.  The changes include all those listed for the various
release candidates (see below) plus the following:

* Several back-end issues are fixed:
  {issue}`388`
  {issue}`389`
  {issue}`390`
  {issue}`392`
* {func}`fwdpy11.TableCollection.fs` added.  See {ref}`tablefs`.
  PR {pr}`387`
  PR {pr}`399`
* Creating populations from `msprime` input improved.
  PR {pr}`395`
* Added {class}`PendingDeprecationWarning` to `fwdpy11.evolve_genomes`.
  PR {pr}`396`

:::{note}

This is the first stable release with support for flexible demographic modeling.
See `softselection` for details as well as `IMexample`.  Currently,
support for different fitness effects in different demes is limited, which
will be addressed in 0.7.0.  However, this version does support adaptation
of quantitative traits to different optima.  See `localadaptation`.

:::

## 0.6.0rc2

Third release candidate of version 0.6.0!

Kind of a big release:

* Fixes a bug in the mechanics of generating offspring metadata.  The bug doesn't
  affect anyone not using custom "genetic value" calculations.  {issue}`371`
* Big reductions in memory requirements for simulations with tree sequence recording.
  PR {pr}`383`
* Better defaults for models with migration.
  PR {pr}`376`
  PR {pr}`375`
  PR {pr}`370`
* Improvements to the C++ back-end of demographic models
  PR {pr}`379`
  PR {pr}`368`
  PR {pr}`367`
  PR {pr}`366`
* Add {class}`fwdpy11.DemographyDebugger`
  PR {pr}`384`
* Add some pre-computed demographic models, see `demographic-models`.
* New examples added:
  `IMexample`
* Many improvements/additions to the test suite and the manual.

## 0.6.0rc1

This is the same as 0.6.0rc0 except that it is based on a master
branch that's been rebased to have the bug fixes from 0.5.5 included.

## 0.6.0rc0

Support for demographic events involving discrete demes.   This is a release
candidate with minimal documentation beyond the examples (see below).

API changes:

* `fwdpy11.Node.population` renamed {attr}`fwdpy11.Node.deme` PR {pr}`340`

This API change won't affect anyone because previous versions didn't support individuals
in different demes.

New features:

* Support for `fwdpy11.DiscreteDemograpy` in simulations with tree sequences.
  PR {pr}`342`
  PR {pr}`346`
  PR {pr}`358`
* Support for different genetic value functions in different demes.
  PR {pr}`357`

Miscellaneous changes:

* Improve how tree sequence nodes are retrieved for "alive" individuals during simulation.
  PR {pr}`344`

New documentation

* Examples of simulations using the `fwdpy11.DiscreteDemograpy` classes.
  PR {pr}`359`
  See `localadaptation` and `migtest`.

Changes to the build system and dependencies:

* Minimum pybind11 version is 2.4.3
* The `-Weffc++` flag is now optional during compilation.

## 0.5.5

This release fixes a rather serious bug.

* Fixes  {issue}`362`
* Fixes  {issue}`363`

The latter is the bad one.  For workflows involving simulate, write
to file, read in and add neutral mutations, that results may now differ.
In practice, we've seen few cases where that has happened (1 in about 10,0000
simulations), but the bug was due to not properly populating a lookup table
of mutation positions after reading the simulation back in from disk.  Thus,
there is the chance that the procedure of putting down neutral mutations
now differs.

## 0.5.4

Bug fix release.

* Fixes  {issue}`350`

## 0.5.3

New features:

* Allow neutral mutations *during* simulations with tree sequences. PR {pr}`328`
* Add C++ back end and Python classes for discrete demographic events. PR {pr}`237`

Miscellaneous changes:

* Links in the manual are now validated via CI. PR {pr}`331`

## 0.5.2

The following bugs are fixed:

* Mutations were not being recycled properly during simulations with tree sequences, resulting in excessive memory consumption. PR {pr}`317`
* Several interface issues with `fwdpy11.MultivariateGSSmo` are fixed. PR {pr}`313`
* Fix a bug that could lead to fixations with tree sequences not "pruning" selected fixations when that behavior is desired. {issue}`287`, fixed in PR {pr}`289`
* A memory safety issue was fixed in the implementation of {attr}`fwdpy11.TreeIterator.samples_below`. PR {pr}`300`.  {issue}`299`

The following new features are added:

* {class}`fwdpy11.BinomialInterval` PR {pr}`322`.
* Allow for preserved samples to be "forgotten" during tree sequence simulations. PR {pr}`306`.

Several performance fixes:

* Extinct genomes are purged at the end of simulations with tree sequences. PR {pr}`319`.
* Improve algorithm to purge extinct variants at the end of a simulation with tree sequences. PR {pr}`318`.
* {func}`fwdpy11.infinite_sites` now returns earlier if possible {issue}`293`.
* Improve performance of mutation counting with ancient samples PR {pr}`289`.

## 0.5.1

This release fixes three bugs:

* `fwdpy11.IndexedEdge` is now exposed to Python. Previously, attempting to access `fwdpy11.TableCollection.input_left` or `fwdpy11.TableCollection.output_right` would give an error because the class contained in these lists wasn't visible. PR {pr}`266`
* {func}`fwdpy11.TreeIterator.roots` now returns the array of roots on the current tree.  Previously, empty arrays were returned. PR {pr}`267`
* Corruption of the samples list using the standalone simplify function. PR {pr}`270`

The following features are new:

* A streamlined API to traverse samples at different time points using {func}`fwdpy11.DiploidPopulation.sample_timepoints`. PR {pr}`279`
* {class}`fwdpy11.TreeIterator` now allows iteration over sites and mutations in the current tree via {func}`fwdpy11.TreeIterator.sites` and {func}`fwdpy11.TreeIterator.mutations`. PR {pr}`275`
* Preorder traversal of nodes in the current tree is possible via {func}`fwdpy11.TreeIterator.nodes`.  Added {func}`fwdpy11.TreeIterator.samples` and {func}`fwdpy11.TreeIterator.samples_below`. PR {pr}`272`

## 0.5.0

This is an intermediate release as we are still working towards supporting more general demographic models.

Major changes include:

* Updating the fwdpp back-end to the pre-release code for fwdpp 0.8.0.  Almost none of these changes are "user facing".
* Add {class}`fwdpy11.SiteTable`, {class}`fwdpy11.Site` and new fields to {class}`fwdpy11.MutationRecord`. PR {pr}`258`  These changes affect the API for some function calls. See {ref}`upgrade_path` for details.

Even though this release changes some of the tree sequence data structures, we are still able to read in files generated by version 0.4.5! (This is actually unit tested.)

Minor changes include:

* Add `fwdpy11.gsl_version`. PR {pr}`256`
* {attr}`fwdpy11.Mutation.g` is converted to the mutation's age when dumping table collections to tskit's format. PR {pr}`257`
* New exception types from fwdpp registered as Python exceptions. PR {pr}`260`
* Several updates to documentation and to continuous integration testing.

## 0.4.5

* {class}`fwdpy11.DataMatrixIterator` now correctly handles nested window coordinates. PR {pr}`244`.

## 0.4.4

* Add {class}`fwdpy11.DataMatrixIterator`. PR {pr}`243`.
* Reduce time needed to execute unit tests of tree sequence functions.

## 0.4.3

* Minor fixes to packaging of source distrubition.
* Add a YCM config file to source repo
* Allow mutation and recombination regions to be empty. PR {pr}`239`.

## 0.4.2

Minor release:

* {class}`fwdpy11.VariantIterator`  may now skip neutral or selected sites during iteration. The behavior is specified
  by parameters passed to the class upon construction.
* Documentation updates

## 0.4.1

Minor release:

* Added position ranges to tree traversal.  PR {pr}`232`.
* Changed default type for range arguments for VariantIterator and data matrix generation. PR {pr}`233`.
* Skipping fixations is now optional in {func}`fwdpy11.data_matrix_from_tables`.
* The C++ back-end for population classes was changed to avoid deleting move constructors. PR {pr}`231`.
* Documentation updates

## 0.4.0

This is a major refactoring:

* The package is now contained in a single namespace, `fwdpy11`.
* The `MlocusPop` concept from previous versions is removed, and {class}`fwdpy11.DiploidPopulation` is the only
  population class now.
* Many Python class names are changed to reflect that there is only one population type now.
* The manual has been rewritten.

The details for this release are best tracked via the cards in [Project 9](https://github.com/molpopgen/fwdpy11/projects/9) on GitHub.

## 0.3.1

Minor bugfix release:

* Preserved nodes are now recorded as samples when table collections are saved to `tskit`
* The fwdpp submodule is updated to include fixes to some debugging code
* Minor updates to the C++ backend of VariantIterator

## 0.3.0

### Deprecations of note

* `fwdpy11.MlocusPop` is *tentatively* deprecated.  The new features described in `geneticmapunit` make
  this class obsolete, but we will await a final verdict pending more testing.

### Bug fixes

* A bug in handling fixations during simulations with tree sequence recording is fixed. This bug is
  GitHub {issue}`200` and the fix is
  PR {pr}`201`.
* Updates to the fwdpp submodule fix a bug in `fwdpy11.ts.infinite_sites`.  Previously, if the genome size
  was not 1.0, then the number of mutations would be off by a factor of the genome size divided by 1.0.  The error was
  due to a bug upstream in fwdpp.
* A bug in how diploid metadata were updated by genetic value types has been fixed.  It is unlikely that this bug
  affected anyone unless they had written custom genetic value calculations where the offspring's genetic value
  depended on the parental metadata. PR {pr}`173`.

### Support for multivariate mutational effects

PR {pr}`164` introduced support for multidimensional mutational effects.
This pull request introduced several changes:

The following new types are added:

* {class}`fwdpy11.MultivariateGaussianEffects`, which is a new "region" type
* `fwdpy11.genetic_values.SlocusPopMultivariateGeneticValueWithMapping`, which is a new ABC for multivariate genetic values
* `fwdpy11.genetic_values.MultivariateGeneticValueToFitnessMap`, which is a new ABC mapping multivariate trait values down to a (single) fitness value.
* `fwdpy11.genetic_values.MultivariateGSS`, which is GSS based on the Euclidean distance from multiple optima
* `fwdpy11.genetic_values.MultivariateGSSmo`, which is the multi-dimensional analog to the existing GSSmo
* `fwdpy11.genetic_values.SlocusMultivariateEffectsStrictAdditive`, which is a new genetic value class for pleiotropic traits.

PR {pr}`175` adds tracking of genetic values during simulation as numpy
arrays via `fwdpy11.Population.genetic_values` and `fwdpy11.Population.ancient_sample_genetic_values`.
Currently, filling these arrays is only supported for simulations with tree sequence recording.

Changes to the C++ back end:

* The API for the C++ class fwdpy11::SlocusPopGeneticValue was slightly changed in order to accommodate the new types.  The old operator() is renamed calculate_gvalue().
* Analogous changes were made to fwdpy11::MlocusPopGeneticValue.

### Dependency changes

* Change minimum GSL version required to 2.3

### Other changes in this release include

It may be helpful to look at the following documentation pages:

* `savingsimstodisk`
* `geneticmapunit`

Detailed changes:

* Add new function to pickle populations while using less memory. PR {pr}`195`,
  PR {pr}`201`
* Improved performance of simulations tracking lots of ancient samples. PR {pr}`194`
* Generalized genetic maps for single-locus simulations.  You can now do much of the "multi-locus" stuff with
  `SlocusPop` now. PR {pr}`189`
* Tree sequence recording now possible for mulit-locus simulations. PR {pr}`185`
* `fwdpy11.ts.count_mutations` added. PR {pr}`183`, PR {pr}`196`, PR {pr}`199`
* Position and key properties added to `fwdpy11.ts.VariantIterator`. PR {pr}`180`
  PR {pr}`181`
* `fwdpy11.ts.TreeIterator` is added, which provides much faster tree traversal. PR {pr}`176`,
  PR {pr}`177`
* `fwdpy11.ts.simplify` no longer retains ancient samples present in the input by default. To do so, explicitly
  label any ancient samples to retain as part of the the samples list passed to the function.
  PR {pr}`169`
* The types {class}`fwdpy11.Region` and {class}`fwdpy11.Sregion` have be re-implemented as C++-based classes, replacing
  the previous pure Python classes.  PR {pr}`163`,
  PR {pr}`174`
* `fwdpy11.model_params.ModelParams.nregions` now defaults to an empty list, which simplifies setup for simulations
  with tree sequences. {commit}`b557c4162cbfdfba6c9126ebec14c7f3f43884eb`.
* When simulating with tree sequences, it is no longer an error to attempt to record ancient samples from the last
  generation of a simulation. PR {pr}`162`

Changes to the C++ back-end include:

* The genetic value types now store a vector of genetic values.  The idea is to generalize the type to handle both uni-
  and multi- variate genetic values. PR {pr}`172`

## Version 0.2.1

This is a point release fixing some minor packaging problems in 0.2.0.

## Version 0.2.0

This release represents major changes to the calclations of genetic values and to how simulations are parameterized.
Please see {ref}`upgrade_path`.

The major feature addition is support for tree sequence recording.  See {ref}`ts-data-types` and {ref}`tsoverview` for details.

### Warning:

This version breaks pickle format compatibility with files generated with version 0.1.4 and earlier.  Sorry, but we had to do it.

### Dependency changes:

* GSL >= 2.2 is now required.
* cmake is now required to build the package.

### Bug fixes:

* Fixed bug in `fwdpy11.util.sort_gamete_keys`.  The function was working on a copy, meaning data were not being
  modified. PR {pr}`93`
* Fix a bug in updating a population's mutation lookup table. This bug was upstream in fwdpp ([fwdpp issue 130](https://github.com/molpopgen/fwdpp/issues/130)).  While definitely a bug, I could never find a case where simulation outputs were adversely affected.  In other words, simulation output remained the same after the fix, due to the rarity of the bug. PR {pr}`98`

### API changes/new features:

* Added support for tree sequence recording.  PR {pr}`142`
* Populations may now be dumped/loaded to/from files. See `fwdpy11.SlocusPop.dump_to_file` and
  `fwdpy11.SlocusPop.load_from_file`.  Analagous functions exist for MlocusPop. PR {pr}`148`
* `fwdpy11.SlocusPop.sample` and `fwdpy11.MlocusPop.sample` now return a `fwdpy11.sampling.DataMatrix`.
  PR {pr}`118`
* `fwdpy11.sampling.DataMatrix` is refactored to match updates to fwdpp.  PR {pr}`139`
* `fwdpy11.sampling.matrix_to_sample` now return a tuple with the neutral and selected data, respectively, as the
  two elements.  PR {pr}`128`
* Diploids have been refactored into two separate classes, {class}`fwdpy11.DiploidGenotype` and
  {class}`fwdpy11.DiploidMetadata`.  Both classes are valid NumPy dtypes.  PR {pr}`108`
* `fwdpy11.model_params.ModelParams` is massively simpilfied. There is now only one class! See `model-params`. PR {pr}`108`
* The design of objects related to calculating genetic values is vastly simplified. PR {pr}`108`
* Populations now contain functions to add mutations, replacing previous functions in fwdpy11.util.  PR {pr}`94`
* `fwdpy11.MlocusPop` now requires that `fwdpy11.MlocusPop.locus_boundaries` be initialized upon
  construction. PR {pr}`96`
* The mutation position lookup table of a population is now a read-only property. PR {pr}`103`
* The mutation position lookup table is now represented as a dict of lists. PR {pr}`121`
* A mutation or fixation can now be rapidy found by its "key".  See `fwdpy11.Population.find_mutation_by_key`
  and `fwdpy11.Population.find_fixation_by_key`.  PR {pr}`106`

### Back-end changes

* The build system now uses cmake.  PR {pr}`151` and {pr}`152`
* Most uses of C's assert macro are replaced with c++ exceptions.  PR {pr}`141`
* The C++ back-end of classes no longer contain any Python objects. PR {pr}`114`
* PR {pr}`108` changes the back-end for representing diploids and for
  calculating genetic values.
* PR {pr}`98` changes the definition of the populaton lookup table, using
  the same model as [fwdpp PR #132](https://github.com/molpopgen/fwdpp/pull/132)
* Refactored class hierarchy for populations. :pr`85`
* Updated to the fwdpp 0.6.x API and cleanup various messes that resulted. PR {pr}`76`, PR {pr}`84`, PR {pr}`90`, PR {pr}`109`, PR {pr}`110`
* The position of extinct variants is set to the max value of a C++ double. PR {pr}`105`
* An entirely new mutation type was introduced on the C++ side.  It is API compatible with the previous type (fwdpp's
  "popgenmut"), but has extra fields for extra flexibility. PR {pr}`77`, PR {pr}`88`
* Replaced `std::bind` with lambda closures for callbacks. PR {pr}`80`
* Fast exposure to raw C++ buffers improved for population objects. PR {pr}`89`
* Refactored long unit tests. PR {pr}`91`
* The GSL error handler is now turned off when fwdpy11 is imported and replaced with a custom handler to propagate GSL errors to C++ exceptions. PR {pr}`140`
* Population mutation position lookup table changed to an unordered multimap. PR {pr}`102`
* When a mutation is fixed or lost, its position is now set to the max value of a C++ double.  This change gets rid of
  some UI oddities when tracking mutations over time. PR {pr}`106` and
  this {commit}`96e8b6e7ca4b257cb8ae5e704f6a36a4b5bfa7bc`.

## Version 0.1.4

### Bug fixes:

* A bug affecting retrieval of multi-locus diploid key data as a buffer for numpy arrays is now fixed. PR {pr}`72`
* `fwdpy11.SingleLocusDiploid.label` is now pickled. PR {pr}`34`

### API changes/new features:

* Population objects have new member functions `sample` and `sample_ind`.  These replace
  `fwdpy11.sampling.sample_separate`, which is now deprecated.  For example, see
  `~fwdpy11.SlocusPop.sample` for more info. (The
  same member functions exist for *all* population objects.) PR {pr}`62`
* Improved support for pickling lower-level types. See the unit test file `tests/test_pickling.py` for examples of directly pickling things like mutations and containers of mutations.  PR {pr}`55`
* `__main__.py` added.  The main use is to help writing python modules based on fwdpy11. See {ref}`developersguide` for details. PR {pr}`54`
* Attributes `popdata` and `popdata_user` added to all population objects. PR {pr}`52`
* `fwdpy11.SingleLocusDiploid.parental_data` added as read-only field. PR {pr}`51`
* `fwdpy11.MlocusPop.locus_boundaries` is now writeable.
* `fwdpy11.sampling.DataMatrix.neutral` and `fwdpy11.sampling.DataMatrix.selected` are now writeable
  buffers. `fwdpy11.sampling.DataMatrix.ndim_neutral` and `fwdpy11.sampling.DataMatrix.ndim_selected` have
  been changed from functions to read-only properties. PR {pr}`45`
* The 'label' field of {class}`fwdpy11.Region` (and {class}`fwdpy11.Sregion`) now populate the label
  field of a mutation. PR {pr}`32` See tests/test_mutation_labels.py for an example.
* Population objects may now be constructed programatically. PR {pr}`36`

### Back-end changes

* The numpy dtype for {class}`fwdpy11.Mutation` has been refactored so that it generates tuples useable to construct object instances. This PR also removes some helper functions in favor of C++11 uniform initialization for these dtypes. PR {pr}`72`
* The documentation building process is greatly streamlined.  PR {pr}`60`
* Object namespaces have been refactored.  The big effect is to streamline the manual. PR {pr}`59`
* Travis CI now tests several Python versions using GCC 6 on Linux. PR {pr}`44`
* `fwdpy11.wright_fisher_qtrait.evolve` has been updated to allow "standard popgen" models of multi-locus
  evolution. This change is a stepping stone to a future global simplification of the API. PR {pr}`42`
* The {class}`fwdpy11.Sregion` now store their callback data differently.  The result is a type that can be
  pickled in Python 3.6. PR {pr}`39`
* Travis builds are now Linux only and test many Python/GCC combos. PR {pr}`38`
* Update to `fwdpp` 0.5.7  PR {pr}`35`
* The method to keep fixations sorted has been updated so that the sorting is by position and fixation time. PR {pr}`33`
* The doctests are now run on Travis. PR {pr}`30`
* Removed all uses of placement new in favor of pybind11::pickle. PR {pr}`26`.
* fwdpy11 are now based on the @property/@foo.setter idiom for safety and code reuse.  PR {pr}`21`

## Version 0.1.3.post1

* Fixed {issue}`23` and {issue}`25` via PR {pr}`24`.

## Version 0.1.3

### Bug fixes:

* {issue}`2` on GitHub fixed. {commit}`562a4d31947d9a7aae31f092ed8c014e94dc56db`

### API changes/new features:

* {class}`fwdpy11.Sregion` may now model distrubitions of effect sizes on scales other than the effect size itself.  A scaling parameter allows the DFE to be functions of N, 2N, 4N, etc. [PR {pr}`16`]
  * Github issues 7, 8, and 9 resolved. All are relatively minor usability tweaks.
* `fwdpy11.util.change_effect_size` added, allowing the "s" and "h" fields of {class}`fwdpy11.Mutation` to be changed. {commit}`ba4841e9407b3d98031801d7eea92b2661871eb2`.
* The attributes of {class}`fwdpy11.Mutation` are now read-only, addressing {issue}`5` on GitHub. {commit}`f376d40788f3d59baa01d1d56b0aa99706560011`
* Trait-to-fitness mapping functions for quantitative trait simulations now take the entire population, rather than just the generation.  This allows us to model things like truncation selection, etc. {commit}`fa37cb8f1763bc7f0e64c8620b6bc1ca350fddb9`

### Back-end changes

* Code base updated to work with [pybind11][pybind11] 2.2.0. [PR {pr}`19`]
* `fwdpy11.model_params` has been refactored, addressing {issue}`4`.  The new code base is more idiomatic w.r.to Python's OO methods. {commit}`1b811c33ab394ae4c64a3c8894984f320b870f22`
* Many of the C++-based types can now be pickled, making model parameter objects easier to serialize.  Most of the
  changes are in {commit}`d0a3602e71a866f7ff9d355d62953ea00c663c5a`.  This mostly addresses {issue}`3`
* Added magic numbers to keep track of compatibility changes to serialization formats.
* __str__ changed to __repr__ for region types {commit}`2df859dd74d3de79d941a1cc21b8712a52bcf9ba`
* fwdpy11.model_params now uses try/except rather than isinstance to check that rates are float-like types. {commit}`37112a60cd8fc74133945e522a47183314bf4085`

## Version 0.1.2

### Bug fixes:

* Fixed bug in setting the number of loci after deserializing a multi-locus population object. {commit}`4e4a547c5b4d30692b62bb4b4a5c22a4cd21d0fa`

### API and back-end changes:

* The C++ data structures are connected to NumPy via Python buffer protocol.  {commit}`48e3925a867c4ec55e1e5bb05457396fb456bc47`
* `fwdpy11.sampling.separate_samples_by_loci` changed to take a list of positions as first argument, and not a population object.

## Version 0.1.1

### Bug fixes:

* Fixed bug in `fwdpy11.sampling.DataMatrix.selected` that returned wrong data in best case scenario and could
  have caused crash in worst case. {commit}`e715fb74472555aa64e1d894563ec218ebba1a97`.
* Fix bug recording fixation times.  If a population was evolved multiple times, fixation times from the later rounds of
  evolution were incorrect. {commit}`9db14d8b3db1c744045e20bfc00ce37e7fb28dfb`
* Fix {issue}`1`, related to fixations in quantitative trait sims. {commit}`6a27386498f056f0c4cc1fc6b8ea12f2b807636c`
* The "label" field of a diploid is now initialized upon constructing a population.

### API and back-end changes:

* Added `fwdpy11.sampling.matrix_to_sample` and `fwdpy11.sampling.separate_samples_by_loci`. {commit}`639c8de999679140fad6a976ff6c1996b25444aa`
* Custom stateless fitness/genetic value calculations may now be implemented with a minimal amount of C++ code.  {commit}`a75166d9ff5471c2d18d66892f9fa01ebec5a667`
* Custom fitness/genetic value calculations now allowed in pure Python, but they are quite slow (for now). {commit}`5549286046ead1181cba684464b3bcb19918321e`
* Stateful trait value models enabled for qtrait sims. {commit}`161dfcef63f3abf28ad56df33b84a92d87d7750f`
* Refactor evolution functions so that stateful fitness models behave as expected.  Enable compiling in a debug mode.
  Fix bug in operator== for diploid type. {commit}`a726c0535a5176aab1df5211fee7bf0aeba5054b`
* fwdpy11.util added, providing `fwdpy11.util.add_mutation`. {commit}`17b92dbe61ee85e2e60211e7dc0ed507a70dbd64`
* Simulations now parameterized using classes in fwdpy11.model_params. {commit}`18e261c8596bf63d2d4e1ef228effb87397b793e` and {commit}`eda7390adb9a98a5d96e6557ba1003488ebac511`
* Added multi-locus simulation of quantitative traits. {commit}`fcad8de9d37bcef5a71ba6d26b4e40e1b67b1993`
* Refactoring of type names. {commit}`632477c7b7592d956149a0cf44e4d26f2a67797e`
* Refactoring internals of single-region fitness/trait value types. {commit}`d55d63631d02fdb2193940475dbcffaa201cf882`
* Allow selected mutations to be retained in fwdpy11.wright_fisher.evolve_regions_sampler_fitness. {commit}`dcc1f2f6555eeada669efef8317f446e3cd0e46a`

**Note:** the refactoring of type names will break scripts based on earlier versions.  Sorry, but things are rapidly changing here.  Please note that you can reassign class and function names in Python, allowing quick hacks to preserve compatibility:

```{code-block} python

import fwdpy11

Spop = fwdpy11.SlocusPop

```

Alternately:

```{code-block} python

from fwdpy11 import SlocusPop as Spop

```

[pybind11]: https://github.com/pybind/pybind11


## 0.19.0a2
