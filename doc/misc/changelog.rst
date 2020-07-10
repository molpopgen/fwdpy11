Changelog
====================================================================================

Major changes are listed below.  Each release likely contains fiddling with back-end code,
updates to latest `fwdpp` version, etc.

0.8.3
****************************************

* :func:`fwdpy11.DiploidPopulation.dump_tables_to_tskit` now populates
  the provenance table. :pr:`542`
* Improve checking migration rates in :class:`fwdpy11.DemographyDebugger`. :pr:`545`
* :class:`fwdpy11.DemographyDebugger` now makes a deep copy of input. :pr:`546`
* The C++ back-ends for Gaussian stabilizing selection classes got streamlined
  without changing the user interface. :pr:`547`
* Manual got overhauled. :pr:`543`
* Snowdrift example (:ref:`here <stateful_fitness>`) reimplemented
  using ``attrs``. :pr:`548`


0.8.2
****************************************

* Fix issue where :class:`fwdpy11.DemographyDebugger` failed to
  catch populations with empty migration matrix rows after
  mass migration. :pr:`539`
* :class:`fwdpy11.DemographyDebugger` is now implemented
  with `attrs`. :pr:`540`.  This change changes a keyword
  argument for this class.  See :ref:`upgrade guide <upgrade_path>`.

0.8.1
****************************************

* Fixed a back-end bug that could have led to corrupt sample lists for simplification. :pr:`536`.
* Made improvements to memory handling of data structures when simulations end. :pr:`537`.
* Added the three-deme model of Jouganous et al. (2017).
  See :func:`fwdpy11.demographic_models.human.jouganous_three_deme`.
  :pr:`534`

0.8.0
****************************************

Stable release. In addition to what is in the previous alpha releases:

* Memory use is substantially reduced due to some improvements
  in `fwdpp`.  :pr:`533` brings in two changes from `fwdpp`,
  :pr:`molpopgen/fwdpp#287` and :pr:`molpopgen/fwdpp#288`.

This release includes some minor `API` changes.
See the :ref:`upgrade guide <upgrade_path>` for details.

0.8.0a1
++++++++++++++++++++++++++++++++++++++++

Second alpha release of 0.8.0:

* Update the ``fwdpp`` submodule. :pr:`529`
* Update documentation related to genetic maps. :pr:`530`

0.8.0a0
++++++++++++++++++++++++++++++++++++++++

This is the first alpha release of 0.8.0.

In addition to what is below, this release contains
a smattering of build system changes, documentation changes,
etc., that are collected under the 0.8.0 milestone on `Github`.

`API`/`UI` changes:

This release brings Python classes that have been reimplemented using `attrs <https://www.attrs.org>`_.  These changes have a lot of benefits:

* A lot of C++ code got removed (yay!) because we can use `attrs` for the pickling
  machinery, `__repr__`, etc..
* We now get *much* nicer `__repr__` for all of the types that get sent into
  instances of :class:`fwdpy11.ModelParams`.

However, these changes required some simplification to the `__init__` methods,
which meant some `API` breakage. See the :ref:`upgrade guide <upgrade_path>`
for details.

This release also removes features deprecated in previous releases. :pr:`482`

Performance improvements:

* Sorting edge tables prior to tree sequence simplification has been replaced 
  by an efficient buffering algorithm. :pr:`526`.

New demographic models:

* The [Tennessen2012]_ model is added via :func:`fwdpy11.demographic_models.human.tennessen`.
  :pr:`479`

Improved behavior:

* Improved warnings about demographic events scheduled to happen
  before the population's current generation. :pr:`495`
* Built-in demographic models now return instances of 
  :class:`fwdpy11.demographic_models.DemographicModelDetails`.
  Such instances can be passed as the `demography` keyword argument
  to initialize :class:`fwdpy11.ModelParams`.
  :pr:`509`.
* The "individual" column of a node table is now populated
  when exporting to a :class:`tskit.TableCollection`. :pr:`488`

Changes to implementation of Python classes

* :class:`fwdpy11.ModelParams` has been reimplemented
  using `attrs <https://www.attrs.org>`_. :pr:`484`, :pr:`486`, :pr:`487`.
* Demographic model types are now implemented using `attrs <https://www.attrs.org>`_ and
  inherit from the C++ back-end class. :pr:`492`
* Region types are now implemented using `attrs <https://www.attrs.org>`_ and
  inherit from the C++ back-end class. :pr:`497`
* Genetic value types are now implemented using `attrs <https://www.attrs.org>`_ and
  inherit from the C++ back-end class. :pr:`504`
* Genetic map unit types are now implemented using `attrs <https://www.attrs.org>`_ and
  inherit from the C++ back-end class. :pr:`506`

C++ back end changes:

* The default C++ language standard is now C++14. :pr:`517`.
* Custom exceptions now have default symbol visibility. :pr:`519`.
* The back-end code for discrete demography got cleaned up. :pr:`521`.
* The `fwdpp` submodule was updated a few times. 
  :pr:`489` :pr:`523` :pr:`525`

0.7.1
++++++++++++++++++++++++++++++++++++++++

Maintenance release and one new feature:

* Allow the first generation of a simulation to be preserved. PR :pr:`470` 
  See :ref:`recapitation`.
* Parameterizing classes like :class:`fwdpy11.GSSmo` is now more Pythonic,
  and some existing ``init`` methods are deprecated in favor of the
  new approach. PR :pr:`461`.

This release include several other improvements to documentation and user interface.
All changes are backwards-compatible, and deprecation warnings are issued when
necessary.  See the 0.7.1 milestone on ``GitHub`` for details.

0.7.0
++++++++++++++++++++++++++++++++++++++++

Major feature release allowing mutations to have different
effect sizes in different demes.

Bugs fixed:

* Temporal samplers now get the correct offspring metadata in simulations
  with tree sequence recording. :issue:`444`

New features:

* Added :class:`fwdpy11.mvDES`, which allows multivariate distributions of effect sizes
  such that mutations have different effect sizes in different demes. See :ref:`mvdes`
  for details. PR :pr:`443` PR :pr:`452`
* :class:`fwdpy11.GeneticValueToFitnessMap` now records whether or not genetic
  values are mapped to fitness or are a trait value via :attr:`fwdpy11.GeneticValueToFitnessMap.maps_to_fitness`
  and :attr:`fwdpy11.GeneticValueToFitnessMap.maps_to_trait_value`.
  PR :pr:`447`

Other changes (see the 0.7.0 milestone on GitHub)

* This release deprecates several features that are no longer sensible given that most
  simulations will use tree sequence recording.  You will see warnings pop up if you
  use these features (or run the unit tests).  These features will be removed
  in 0.8.0.
* Many back-end changes to the C++ code simplify things in various places.

0.6.4
++++++++++++++++++++++++++++++++++++++++

Fixes a bug where the timing of updates to stateful genetic values
was off by one generation:

*  :issue:`437`

0.6.3
++++++++++++++++++++++++++++++++++++++++

Maintenance release.

This release fixes three bugs. The first two are related to internal
details of book-keeping various data structures:

*  :issue:`420`
*  :issue:`422`
*  :issue:`432`

Other changes:

* :mod:`sparse` is added to ``install_requires`` in ``setup.py``.  :issue:`421`
* :class:`fwdpy11.TableCollection`'s validation of genome lengths is improved. PR :pr:`428` 
* The C++ base class for a population is now a concrete class rather than a template alias.  This change enables forward declarations in header files. PR :pr:`427` 

0.6.2
++++++++++++++++++++++++++++++++++++++++

This release changes the migration code to model juvenile migration.
These changes simplify the back end and give the same results (in
distribution).  The relevant PRs are:

* PR :pr:`416` 
* PR :pr:`417` 

0.6.1
++++++++++++++++++++++++++++++++++++++++

This is a maintenance release that clears up a few issues:

*  :issue:`246`
*  :issue:`280`
*  :issue:`339`
*  :issue:`365`
*  :issue:`386`
*  :issue:`397`

The following features are added:

* :attr:`fwdpy11.DataMatrix.neutral_matrix`
* :attr:`fwdpy11.DataMatrix.selected_matrix`
* :func:`fwdpy11.DataMatrix.merge`

0.6.0
++++++++++++++++++++++++++++++++++++++++

This is a major feature release.  The changes include all those listed for the various 
release candidates (see below) plus the following:

* Several back-end issues are fixed: 
  :issue:`388`
  :issue:`389`
  :issue:`390`
  :issue:`392`
* :func:`fwdpy11.TableCollection.fs` added.  See :ref:`tablefs`.
  PR :pr:`387` 
  PR :pr:`399` 
* Creating populations from :mod:`msprime` input improved.
  PR :pr:`395` 
* Added :class:`PendingDeprecationWarning` to ``fwdpy11.evolve_genomes``.
  PR :pr:`396` 

.. note::

    This is the first stable release with support for flexible demographic modeling.
    See :ref:`softselection` for details as well as :ref:`IMexample`.  Currently,
    support for different fitness effects in different demes is limited, which
    will be addressed in 0.7.0.  However, this version does support adaptation
    of quantitative traits to different optima.  See :ref:`localadaptation`.


0.6.0rc2 
++++++++++++++++

Third release candidate of version 0.6.0!

Kind of a big release:

* Fixes a bug in the mechanics of generating offspring metadata.  The bug doesn't
  affect anyone not using custom "genetic value" calculations.  :issue:`371`
* Big reductions in memory requirements for simulations with tree sequence recording.
  PR :pr:`383` 
* Better defaults for models with migration.
  PR :pr:`376` 
  PR :pr:`375` 
  PR :pr:`370` 
* Improvements to the C++ back-end of demographic models
  PR :pr:`379` 
  PR :pr:`368` 
  PR :pr:`367` 
  PR :pr:`366` 
* Add :class:`fwdpy11.DemographyDebugger`
  PR :pr:`384` 
* Add some pre-computed demographic models, see :ref:`demographic_models`.
* New examples added:
  :ref:`IMexample`
* Many improvements/additions to the test suite and the manual.
  
  
0.6.0rc1
++++++++++++++++

This is the same as 0.6.0rc0 except that it is based on a master
branch that's been rebased to have the bug fixes from 0.5.5 included.

0.6.0rc0
++++++++++++++++

Support for demographic events involving discrete demes.   This is a release 
candidate with minimal documentation beyond the examples (see below).

API changes:

* ``fwdpy11.Node.population`` renamed :attr:`fwdpy11.Node.deme` PR :pr:`340`

This API change won't affect anyone because previous versions didn't support individuals
in different demes.

New features:

* Support for :class:`fwdpy11.DiscreteDemography` in simulations with tree sequences.
  PR :pr:`342` 
  PR :pr:`346` 
  PR :pr:`358` 

* Support for different genetic value functions in different demes. 
  PR :pr:`357` 

Miscellaneous changes:

* Improve how tree sequence nodes are retrieved for "alive" individuals during simulation.
  PR :pr:`344` 
   
New documentation

* Examples of simulations using the :class:`fwdpy11.DiscreteDemography` classes.
  PR :pr:`359` 
  See :ref:`localadaptation` and :ref:`migtest`.

Changes to the build system and dependencies:

* Minimum pybind11 version is 2.4.3
* The ``-Weffc++`` flag is now optional during compilation.

0.5.5
++++++++++++++++

This release fixes a rather serious bug.

* Fixes  :issue:`362`
* Fixes  :issue:`363`

The latter is the bad one.  For workflows involving simulate, write
to file, read in and add neutral mutations, that results may now differ.
In practice, we've seen few cases where that has happened (1 in about 10,0000
simulations), but the bug was due to not properly populating a lookup table
of mutation positions after reading the simulation back in from disk.  Thus,
there is the chance that the procedure of putting down neutral mutations
now differs.

0.5.4
++++++++++++++++

Bug fix release.

* Fixes  :issue:`350`

0.5.3
++++++++++++++++

New features:

* Allow neutral mutations *during* simulations with tree sequences. PR :pr:`328`
* Add C++ back end and Python classes for discrete demographic events. PR :pr:`237` 

Miscellaneous changes:

* Links in the manual are now validated via CI. PR :pr:`331` 

0.5.2
++++++++++++++++

The following bugs are fixed:

* Mutations were not being recycled properly during simulations with tree sequences, resulting in excessive memory consumption. PR :pr:`317`
* Several interface issues with :class:`fwdpy11.MultivariateGSSmo` are fixed. PR :pr:`313`
* Fix a bug that could lead to fixations with tree sequences not "pruning" selected fixations when that behavior is desired. :issue:`287`, fixed in PR :pr:`289`
* A memory safety issue was fixed in the implementation of :attr:`fwdpy11.TreeIterator.samples_below`. PR :pr:`300`.  :issue:`299`

The following new features are added:

* :class:`fwdpy11.BinomialInterval` PR :pr:`322`.
* Allow for preserved samples to be "forgotten" during tree sequence simulations. PR :pr:`306`. See :ref:`tstimeseries`

Several performance fixes:

* Extinct genomes are purged at the end of simulations with tree sequences. PR :pr:`319`.
* Improve algorithm to purge extinct variants at the end of a simulation with tree sequences. PR :pr:`318`.
* :func:`fwdpy11.infinite_sites` now returns earlier if possible :issue:`293`.
* Improve performance of mutation counting with ancient samples PR :pr:`289`.


0.5.1
++++++++++++++++

This release fixes three bugs:

* ``fwdpy11.IndexedEdge`` is now exposed to Python. Previously, attempting to access `fwdpy11.TableCollection.input_left` or `fwdpy11.TableCollection.output_right` would give an error because the class contained in these lists wasn't visible. PR :pr:`266`
* :func:`fwdpy11.TreeIterator.roots` now returns the array of roots on the current tree.  Previously, empty arrays were returned. PR :pr:`267`
* Corruption of the samples list using the standalone simplify function. PR :pr:`270`

The following features are new:

* A streamlined API to traverse samples at different time points using :func:`fwdpy11.DiploidPopulation.sample_timepoints`. PR :pr:`279`
* :class:`fwdpy11.TreeIterator` now allows iteration over sites and mutations in the current tree via :func:`fwdpy11.TreeIterator.sites` and :func:`fwdpy11.TreeIterator.mutations`. PR :pr:`275`
* Preorder traversal of nodes in the current tree is possible via :func:`fwdpy11.TreeIterator.nodes`.  Added :func:`fwdpy11.TreeIterator.samples` and :func:`fwdpy11.TreeIterator.samples_below`. PR :pr:`272`

0.5.0
+++++++++++

This is an intermediate release as we are still working towards supporting more general demographic models.

Major changes include:

* Updating the fwdpp back-end to the pre-release code for fwdpp 0.8.0.  Almost none of these changes are "user facing".
* Add :class:`fwdpy11.SiteTable`, :class:`fwdpy11.Site` and new fields to :class:`fwdpy11.MutationRecord`. PR :pr:`258`  These changes affect the API for some function calls. See :ref:`upgrade_path` for details.

Even though this release changes some of the tree sequence data structures, we are still able to read in files generated by version 0.4.5! (This is actually unit tested.)

Minor changes include:

* Add `fwdpy11.gsl_version`. PR :pr:`256`
* :attr:`fwdpy11.Mutation.g` is converted to the mutation's age when dumping table collections to tskit's format. PR :pr:`257`
* New exception types from fwdpp registered as Python exceptions. PR :pr:`260`
* Several updates to documentation and to continuous integration testing.

0.4.5
+++++++++++

* :class:`fwdpy11.DataMatrixIterator` now correctly handles nested window coordinates. PR :pr:`244`.


0.4.4
+++++++++++

* Add :class:`fwdpy11.DataMatrixIterator`. PR :pr:`243`.
* Reduce time needed to execute unit tests of tree sequence functions.

0.4.3
++++++++++++++++++++++++++++++++++

* Minor fixes to packaging of source distrubition.
* Add a YCM config file to source repo
* Allow mutation and recombination regions to be empty. PR :pr:`239`.

0.4.2
++++++++++++++++++++++++++++++++++

Minor release:

* :class:`fwdpy11.VariantIterator`  may now skip neutral or selected sites during iteration. The behavior is specified
  by parameters passed to the class upon construction.
* Documentation updates

0.4.1
++++++++++++++++++++++++++++++++++

Minor release:

* Added position ranges to tree traversal.  PR :pr:`232`.
* Changed default type for range arguments for VariantIterator and data matrix generation. PR :pr:`233`.
* Skipping fixations is now optional in :func:`fwdpy11.data_matrix_from_tables`.
* The C++ back-end for population classes was changed to avoid deleting move constructors. PR :pr:`231`.
* Documentation updates

0.4.0
++++++++++++++++++++++++++++++++++

This is a major refactoring:

* The package is now contained in a single namespace, `fwdpy11`.
* The `MlocusPop` concept from previous versions is removed, and :class:`fwdpy11.DiploidPopulation` is the only
  population class now.
* Many Python class names are changed to reflect that there is only one population type now.
* The manual has been rewritten.

The details for this release are best tracked via the cards in `Project 9 <https://github.com/molpopgen/fwdpy11/projects/9>`_ on GitHub.


0.3.1
++++++++++++++++++++++++++

Minor bugfix release:

* Preserved nodes are now recorded as samples when table collections are saved to `tskit`
* The fwdpp submodule is updated to include fixes to some debugging code
* Minor updates to the C++ backend of VariantIterator

0.3.0
++++++++++++++++++++++++++

Deprecations of note
-------------------------------------------------------------

* `fwdpy11.MlocusPop` is *tentatively* deprecated.  The new features described in :ref:`geneticmapunit` make
  this class obsolete, but we will await a final verdict pending more testing.

Bug fixes
-------------------------------------------------------------

* A bug in handling fixations during simulations with tree sequence recording is fixed. This bug is 
  GitHub :issue:`200` and the fix is
  PR :pr:`201`.
* Updates to the fwdpp submodule fix a bug in :func:`fwdpy11.ts.infinite_sites`.  Previously, if the genome size 
  was not 1.0, then the number of mutations would be off by a factor of the genome size divided by 1.0.  The error was
  due to a bug upstream in fwdpp.
* A bug in how diploid metadata were updated by genetic value types has been fixed.  It is unlikely that this bug
  affected anyone unless they had written custom genetic value calculations where the offspring's genetic value 
  depended on the parental metadata. PR :pr:`173`. 

Support for multivariate mutational effects
-------------------------------------------------------------

PR :pr:`164` introduced support for multidimensional mutational effects.
This pull request introduced several changes: 

The following new types are added:

* :class:`fwdpy11.MultivariateGaussianEffects`, which is a new "region" type
* :class:`fwdpy11.genetic_values.SlocusPopMultivariateGeneticValueWithMapping`, which is a new ABC for multivariate genetic values
* :class:`fwdpy11.genetic_values.MultivariateGeneticValueToFitnessMap`, which is a new ABC mapping multivariate trait values down to a (single) fitness value.
* :class:`fwdpy11.genetic_values.MultivariateGSS`, which is GSS based on the Euclidean distance from multiple optima
* :class:`fwdpy11.genetic_values.MultivariateGSSmo`, which is the multi-dimensional analog to the existing GSSmo
* :class:`fwdpy11.genetic_values.SlocusMultivariateEffectsStrictAdditive`, which is a new genetic value class for pleiotropic traits.

PR :pr:`175` adds tracking of genetic values during simulation as numpy
arrays via :attr:`fwdpy11.Population.genetic_values` and :attr:`fwdpy11.Population.ancient_sample_genetic_values`.
Currently, filling these arrays is only supported for simulations with tree sequence recording.

Changes to the C++ back end:

* The API for the C++ class fwdpy11::SlocusPopGeneticValue was slightly changed in order to accommodate the new types.  The old operator() is renamed calculate_gvalue().
* Analogous changes were made to fwdpy11::MlocusPopGeneticValue.


Dependency changes
-------------------------------------------------------------

* Change minimum GSL version required to 2.3

Other changes in this release include
-------------------------------------------------------------

It may be helpful to look at the following documentation pages:

* :ref:`savingsimstodisk`
* :ref:`geneticmapunit`

Detailed changes:

* Add new function to pickle populations while using less memory. PR :pr:`195`,
  PR :pr:`201`
* Improved performance of simulations tracking lots of ancient samples. PR :pr:`194`
* Generalized genetic maps for single-locus simulations.  You can now do much of the "multi-locus" stuff with
  `SlocusPop` now. PR :pr:`189`
* Tree sequence recording now possible for mulit-locus simulations. PR :pr:`185`
* :func:`fwdpy11.ts.count_mutations` added. PR :pr:`183`, PR :pr:`196`, PR :pr:`199`
* Position and key properties added to :class:`fwdpy11.ts.VariantIterator`. PR :pr:`180`
  PR :pr:`181`
* :class:`fwdpy11.ts.TreeIterator` is added, which provides much faster tree traversal. PR :pr:`176`,
  PR :pr:`177`
* :func:`fwdpy11.ts.simplify` no longer retains ancient samples present in the input by default. To do so, explicitly
  label any ancient samples to retain as part of the the samples list passed to the function.
  PR :pr:`169`
* The types :class:`fwdpy11.Region` and :class:`fwdpy11.Sregion` have be re-implemented as C++-based classes, replacing 
  the previous pure Python classes.  PR :pr:`163`,
  PR :pr:`174`
* :attr:`fwdpy11.model_params.ModelParams.nregions` now defaults to an empty list, which simplifies setup for simulations
  with tree sequences. :commit:`b557c4162cbfdfba6c9126ebec14c7f3f43884eb`. 
* When simulating with tree sequences, it is no longer an error to attempt to record ancient samples from the last
  generation of a simulation. PR :pr:`162`

Changes to the C++ back-end include:

* The genetic value types now store a vector of genetic values.  The idea is to generalize the type to handle both uni-
  and multi- variate genetic values. PR :pr:`172`

Version 0.2.1
++++++++++++++++++++++++++

This is a point release fixing some minor packaging problems in 0.2.0.

Version 0.2.0
++++++++++++++++++++++++++

This release represents major changes to the calclations of genetic values and to how simulations are parameterized.
Please see :ref:`upgrade_path`, :ref:`genetic_values_types`, and :ref:`model_params` for details.

The major feature addition is support for tree sequence recording.  See :ref:`ts_data_types` and :ref:`ts` for details.

Warning:
--------------------------

This version breaks pickle format compatibility with files generated with version 0.1.4 and earlier.  Sorry, but we had to do it.

Dependency changes:
--------------------------

* GSL >= 2.2 is now required.
* cmake is now required to build the package.

Bug fixes:
--------------------------

* Fixed bug in :func:`fwdpy11.util.sort_gamete_keys`.  The function was working on a copy, meaning data were not being
  modified. PR :pr:`93`
* Fix a bug in updating a population's mutation lookup table. This bug was upstream in fwdpp (`fwdpp issue 130 <https://github.com/molpopgen/fwdpp/issues/130>`_).  While definitely a bug, I could never find a case where simulation outputs were adversely affected.  In other words, simulation output remained the same after the fix, due to the rarity of the bug. PR :pr:`98`


API changes/new features:
----------------------------------------------------

* Added support for tree sequence recording.  PR :pr:`142`
* Populations may now be dumped/loaded to/from files. See :func:`fwdpy11.SlocusPop.dump_to_file` and
  :func:`fwdpy11.SlocusPop.load_from_file`.  Analagous functions exist for MlocusPop. PR :pr:`148`
* :func:`fwdpy11.SlocusPop.sample` and :func:`fwdpy11.MlocusPop.sample` now return a :class:`fwdpy11.sampling.DataMatrix`.
  PR :pr:`118`
* :class:`fwdpy11.sampling.DataMatrix` is refactored to match updates to fwdpp.  PR :pr:`139`
* :func:`fwdpy11.sampling.matrix_to_sample` now return a tuple with the neutral and selected data, respectively, as the
  two elements.  PR :pr:`128`
* Diploids have been refactored into two separate classes, :class:`fwdpy11.DiploidGenotype` and
  :class:`fwdpy11.DiploidMetadata`.  Both classes are valid NumPy dtypes.  See :ref:`processingpopsNP`. PR :pr:`108`
* :class:`fwdpy11.model_params.ModelParams` is massively simpilfied. There is now only one class! See :ref:`model_params`. PR :pr:`108`
* The design of objects related to calculating genetic values is vastly simplified.  See :ref:`genetic_values_types`. PR :pr:`108`
* Populations now contain functions to add mutations, replacing previous functions in fwdpy11.util.  PR :pr:`94`
* :class:`fwdpy11.MlocusPop` now requires that :attr:`fwdpy11.MlocusPop.locus_boundaries` be initialized upon
  construction. PR :pr:`96`
* The mutation position lookup table of a population is now a read-only property. See :ref:`mpos`. PR :pr:`103`
* The mutation position lookup table is now represented as a dict of lists. PR :pr:`121`
* A mutation or fixation can now be rapidy found by its "key".  See :func:`fwdpy11.Population.find_mutation_by_key`
  and :func:`fwdpy11.Population.find_fixation_by_key`.  PR :pr:`106`

Back-end changes
------------------------

* The build system now uses cmake.  PR :pr:`151` and :pr:`152`
* Most uses of C's assert macro are replaced with c++ exceptions.  PR :pr:`141`
* The C++ back-end of classes no longer contain any Python objects. PR :pr:`114`
* PR :pr:`108` changes the back-end for representing diploids and for
  calculating genetic values.
* PR :pr:`98` changes the definition of the populaton lookup table, using
  the same model as `fwdpp PR #132 <https://github.com/molpopgen/fwdpp/pull/132>`_
* Refactored class hierarchy for populations. :pr`85`
* Updated to the fwdpp 0.6.x API and cleanup various messes that resulted. PR :pr:`76`, PR :pr:`84`, PR :pr:`90`, PR :pr:`109`, PR :pr:`110`
* The position of extinct variants is set to the max value of a C++ double. PR :pr:`105`
* An entirely new mutation type was introduced on the C++ side.  It is API compatible with the previous type (fwdpp's
  "popgenmut"), but has extra fields for extra flexibility. PR :pr:`77`, PR :pr:`88`
* Replaced `std::bind` with lambda closures for callbacks. PR :pr:`80`
* Fast exposure to raw C++ buffers improved for population objects. PR :pr:`89`
* Refactored long unit tests. PR :pr:`91`
* The GSL error handler is now turned off when fwdpy11 is imported and replaced with a custom handler to propagate GSL errors to C++ exceptions. PR :pr:`140`
* Population mutation position lookup table changed to an unordered multimap. PR :pr:`102`
* When a mutation is fixed or lost, its position is now set to the max value of a C++ double.  This change gets rid of
  some UI oddities when tracking mutations over time. PR :pr:`106` and
  this :commit:`96e8b6e7ca4b257cb8ae5e704f6a36a4b5bfa7bc`.

Version 0.1.4
++++++++++++++++++++++++++

Bug fixes:
--------------------------

* A bug affecting retrieval of multi-locus diploid key data as a buffer for numpy arrays is now fixed. PR :pr:`72`
* :attr:`fwdpy11.SingleLocusDiploid.label` is now pickled. PR :pr:`34`
    
API changes/new features:
----------------------------------------------------

* Population objects have new member functions ``sample`` and ``sample_ind``.  These replace
  :func:`fwdpy11.sampling.sample_separate`, which is now deprecated.  For example, see
  :func:`~fwdpy11.SlocusPop.sample` for more info. (The
  same member functions exist for *all* population objects.) PR :pr:`62`
* Improved support for pickling lower-level types. See the unit test file `tests/test_pickling.py` for examples of directly pickling things like mutations and containers of mutations.  PR :pr:`55`
* `__main__.py` added.  The main use is to help writing python modules based on fwdpy11. See :ref:`developers` for details. PR :pr:`54`
* Attributes `popdata` and `popdata_user` added to all population objects. PR :pr:`52`
* :attr:`fwdpy11.SingleLocusDiploid.parental_data` added as read-only field. PR :pr:`51`
* :attr:`fwdpy11.MlocusPop.locus_boundaries` is now writeable.
* :attr:`fwdpy11.sampling.DataMatrix.neutral` and :attr:`fwdpy11.sampling.DataMatrix.selected` are now writeable
  buffers. :attr:`fwdpy11.sampling.DataMatrix.ndim_neutral` and :attr:`fwdpy11.sampling.DataMatrix.ndim_selected` have
  been changed from functions to read-only properties. PR :pr:`45`
* The 'label' field of :class:`fwdpy11.Region` (and :class:`fwdpy11.Sregion`) now populate the label
  field of a mutation. PR :pr:`32` See tests/test_mutation_labels.py for an example.
* Population objects may now be constructed programatically. See :ref:`popobjects`.   PR :pr:`36` 

Back-end changes
------------------------

* The numpy dtype for :class:`fwdpy11.Mutation` has been refactored so that it generates tuples useable to construct object instances. This PR also removes some helper functions in favor of C++11 uniform initialization for these dtypes. PR :pr:`72`
* The documentation building process is greatly streamlined.  PR :pr:`60`
* Object namespaces have been refactored.  The big effect is to streamline the manual. PR :pr:`59`
* Travis CI now tests several Python versions using GCC 6 on Linux. PR :pr:`44`
* :func:`fwdpy11.wright_fisher_qtrait.evolve` has been updated to allow "standard popgen" models of multi-locus
  evolution. This change is a stepping stone to a future global simplification of the API. PR :pr:`42`
* The :class:`fwdpy11.Sregion` now store their callback data differently.  The result is a type that can be
  pickled in Python 3.6. PR :pr:`39` 
* Travis builds are now Linux only and test many Python/GCC combos. PR :pr:`38`
* Update to fwdpp_ 0.5.7  PR :pr:`35`
* The method to keep fixations sorted has been updated so that the sorting is by position and fixation time. PR :pr:`33`
* The doctests are now run on Travis. PR :pr:`30`
* Removed all uses of placement new in favor of pybind11::pickle. PR :pr:`26`.
* fwdpy11 are now based on the @property/@foo.setter idiom for safety and code reuse.  PR :pr:`21`

Version 0.1.3.post1
++++++++++++++++++++++++++

* Fixed :issue:`23` and :issue:`25` via PR :pr:`24`.

Version 0.1.3
++++++++++++++++++++++++++

Bug fixes:
------------------------

* :issue:`2` on GitHub fixed. :commit:`562a4d31947d9a7aae31f092ed8c014e94dc56db`

API changes/new features:
------------------------------------------------

* :class:`fwdpy11.Sregion` may now model distrubitions of effect sizes on scales other than the effect size itself.  A scaling parameter allows the DFE to be functions of N, 2N, 4N, etc. [PR :pr:`16`]
  * Github issues 7, 8, and 9 resolved. All are relatively minor usability tweaks.
* :func:`fwdpy11.util.change_effect_size` added, allowing the "s" and "h" fields of :class:`fwdpy11.Mutation` to be changed. :commit:`ba4841e9407b3d98031801d7eea92b2661871eb2`.
* The attributes of :class:`fwdpy11.Mutation` are now read-only, addressing :issue:`5` on GitHub. :commit:`f376d40788f3d59baa01d1d56b0aa99706560011`
* Trait-to-fitness mapping functions for quantitative trait simulations now take the entire population, rather than just the generation.  This allows us to model things like truncation selection, etc. :commit:`fa37cb8f1763bc7f0e64c8620b6bc1ca350fddb9`

Back-end changes
------------------------

* Code base updated to work with pybind11_ 2.2.0. [PR :pr:`19`] 
* :mod:`fwdpy11.model_params` has been refactored, addressing :issue:`4`.  The new code base is more idiomatic w.r.to Python's OO methods. :commit:`1b811c33ab394ae4c64a3c8894984f320b870f22`
* Many of the C++-based types can now be pickled, making model parameter objects easier to serialize.  Most of the
  changes are in :commit:`d0a3602e71a866f7ff9d355d62953ea00c663c5a`.  This mostly addresses :issue:`3`
* Added magic numbers to keep track of compatibility changes to serialization formats.
* __str__ changed to __repr__ for region types :commit:`2df859dd74d3de79d941a1cc21b8712a52bcf9ba`
* fwdpy11.model_params now uses try/except rather than isinstance to check that rates are float-like types. :commit:`37112a60cd8fc74133945e522a47183314bf4085`

Version 0.1.2
++++++++++++++++++++++++++

Bug fixes:
---------------------
* Fixed bug in setting the number of loci after deserializing a multi-locus population object. :commit:`4e4a547c5b4d30692b62bb4b4a5c22a4cd21d0fa`

API and back-end changes:
------------------------------------------
* The C++ data structures are connected to NumPy via Python buffer protocol.  See :ref:`processingpopsNP`.  :commit:`48e3925a867c4ec55e1e5bb05457396fb456bc47`
* :func:`fwdpy11.sampling.separate_samples_by_loci` changed to take a list of positions as first argument, and not a population object.

Version 0.1.1
++++++++++++++++++++++++++

Bug fixes:
---------------------
* Fixed bug in :func:`fwdpy11.sampling.DataMatrix.selected` that returned wrong data in best case scenario and could
  have caused crash in worst case. :commit:`e715fb74472555aa64e1d894563ec218ebba1a97`.
* Fix bug recording fixation times.  If a population was evolved multiple times, fixation times from the later rounds of
  evolution were incorrect. :commit:`9db14d8b3db1c744045e20bfc00ce37e7fb28dfb`
* Fix :issue:`1`, related to fixations in quantitative trait sims. :commit:`6a27386498f056f0c4cc1fc6b8ea12f2b807636c`
* The "label" field of a diploid is now initialized upon constructing a population.

API and back-end changes:
------------------------------------------
* Added :func:`fwdpy11.sampling.matrix_to_sample` and :func:`fwdpy11.sampling.separate_samples_by_loci`. :commit:`639c8de999679140fad6a976ff6c1996b25444aa`
* Custom stateless fitness/genetic value calculations may now be implemented with a minimal amount of C++ code. See
  :ref:`customgvalues`. :commit:`a75166d9ff5471c2d18d66892f9fa01ebec5a667`
* Custom fitness/genetic value calculations now allowed in pure Python, but they are quite slow (for now). See 
  :ref:`customgvalues`. :commit:`5549286046ead1181cba684464b3bcb19918321e`
* Stateful trait value models enabled for qtrait sims. :commit:`161dfcef63f3abf28ad56df33b84a92d87d7750f`
* Refactor evolution functions so that stateful fitness models behave as expected.  Enable compiling in a debug mode.
  Fix bug in operator== for diploid type. :commit:`a726c0535a5176aab1df5211fee7bf0aeba5054b`
* fwdpy11.util added, providing :func:`fwdpy11.util.add_mutation`. :commit:`17b92dbe61ee85e2e60211e7dc0ed507a70dbd64`
* Simulations now parameterized using classes in fwdpy11.model_params. :commit:`18e261c8596bf63d2d4e1ef228effb87397b793e` and :commit:`eda7390adb9a98a5d96e6557ba1003488ebac511`
* Added multi-locus simulation of quantitative traits. :commit:`fcad8de9d37bcef5a71ba6d26b4e40e1b67b1993`
* Refactoring of type names. :commit:`632477c7b7592d956149a0cf44e4d26f2a67797e`
* Refactoring internals of single-region fitness/trait value types. :commit:`d55d63631d02fdb2193940475dbcffaa201cf882`
* Allow selected mutations to be retained in fwdpy11.wright_fisher.evolve_regions_sampler_fitness. :commit:`dcc1f2f6555eeada669efef8317f446e3cd0e46a`

**Note:** the refactoring of type names will break scripts based on earlier versions.  Sorry, but things are rapidly changing here.  Please note that you can reassign class and function names in Python, allowing quick hacks to preserve compatibility:

.. code-block:: python

    import fwdpy11

    Spop = fwdpy11.SlocusPop

Alternately:

.. code-block:: python

    from fwdpy11 import SlocusPop as Spop

.. _pybind11: https://github.com/pybind/pybind11

