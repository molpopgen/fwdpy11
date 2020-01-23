Changelog
====================================================================================

Major changes are listed below.  Each release likely contains fiddling with back-end code, updates to latest fwdpp
version, etc.

0.6.0rc0
++++++++++++++++

Support for demographic events involving discrete demes.   This is a release 
candidate with minimal documentation beyond the examples (see below).

API changes:

* ``fwdpy11.Node.population`` renamed :attr:`fwdpy11.Node.deme` `PR340 <https://github.com/molpopgen/fwdpy11/pull/340>`_

This API change won't affect anyone because previous versions didn't support individuals
in different demes.

New features:

* Support for :class:`fwdpy11.DiscreteDemography` in simulations with tree sequences.
  `PR342 <https://github.com/molpopgen/fwdpy11/pull/342>`_ 
  `PR346 <https://github.com/molpopgen/fwdpy11/pull/346>`_ 
  `PR358 <https://github.com/molpopgen/fwdpy11/pull/358>`_ 

* Support for different genetic value functions in different demes. 
  `PR357 <https://github.com/molpopgen/fwdpy11/pull/357>`_ 

Miscellaneous changes:

* Improve how tree sequence nodes are retrieved for "alive" individuals during simulation.
  `PR344 <https://github.com/molpopgen/fwdpy11/pull/344>`_ 
   
New documentation

* Examples of simulations using the :class:`fwdpy11.DiscreteDemography` classes.
  `PR359 <https://github.com/molpopgen/fwdpy11/pull/359>`_ 
  See :ref:`localadaptation` and :ref:`migtest`.

Changes to the build system and dependencies:

* Minimum pybind11 version is 2.4.3
* The ``-Weffc++`` flag is now optional during compilation.


0.5.4
++++++++++++++++

Bug fix release.

* Fixes  `PR350 <https://github.com/molpopgen/fwdpy11/issues/350>`_

0.5.3
++++++++++++++++

New features:

* Allow neutral mutations *during* simulations with tree sequences. `PR328 <https://github.com/molpopgen/fwdpy11/pull/328>`_
* Add C++ back end and Python classes for discrete demographic events. `PR237 <https://github.com/molpopgen/fwdpy11/pull/237>`_ 

Miscellaneous changes:

* Links in the manual are now validated via CI. `PR331 <https://github.com/molpopgen/fwdpy11/pull/331>`_ 

0.5.2
++++++++++++++++

The following bugs are fixed:

* Mutations were not being recycled properly during simulations with tree sequences, resulting in excessive memory consumption. `PR317 <https://github.com/molpopgen/fwdpy11/pull/317>`_
* Several interface issues with :class:`fwdpy11.MultivariateGSSmo` are fixed. `PR313 <https://github.com/molpopgen/fwdpy11/pull/313>`_
* Fix a bug that could lead to fixations with tree sequences not "pruning" selected fixations when that behavior is desired. `Issue287 <https://github.com/molpopgen/fwdpy11/issues/287>`_, fixed in `PR289 <https://github.com/molpopgen/fwdpy11/pull/289>`_
* A memory safety issue was fixed in the implementation of :attr:`fwdpy11.TreeIterator.samples_below`. `PR300 <https://github.com/molpopgen/fwdpy11/pull/300>`_.  `Issue299 <https://github.com/molpopgen/fwdpy11/issues/299>`_

The following new features are added:

* :class:`fwdpy11.BinomialInterval` `PR322 <https://github.com/molpopgen/fwdpy11/pull/322>`_.
* Allow for preserved samples to be "forgotten" during tree sequence simulations. `PR306 <https://github.com/molpopgen/fwdpy11/pull/306>`_. See :ref:`tstimeseries`

Several performance fixes:

* Extinct genomes are purged at the end of simulations with tree sequences. `PR319 <https://github.com/molpopgen/fwdpy11/pull/319>`_.
* Improve algorithm to purge extinct variants at the end of a simulation with tree sequences. `PR318 <https://github.com/molpopgen/fwdpy11/pull/318>`_.
* :func:`fwdpy11.infinite_sites` now returns earlier if possible `Issue293 <https://github.com/molpopgen/fwdpy11/issues/293>`_.
* Improve performance of mutation counting with ancient samples `PR289 <https://github.com/molpopgen/fwdpy11/pull/289>`_.


0.5.1
++++++++++++++++

This release fixes three bugs:

* :class:`fwdpy11.IndexedEdge` is now exposed to Python. Previously, attempting to access :attr:`fwdpy11.TableCollection.input_left` or :attr:`fwdpy11.TableCollection.output_right` would give an error because the class contained in these lists wasn't visible. `PR 266 <https://github.com/molpopgen/fwdpy11/pull/266>`_
* :func:`fwdpy11.TreeIterator.roots` now returns the array of roots on the current tree.  Previously, empty arrays were returned. `PR 267 <https://github.com/molpopgen/fwdpy11/pull/267>`_
* Corruption of the samples list using the standalone simplify function. `PR 270 <https://github.com/molpopgen/fwdpy11/pull/270>`_

The following features are new:

* A streamlined API to traverse samples at different time points using :func:`fwdpy11.DiploidPopulation.sample_timepoints`. `PR 279 <https://github.com/molpopgen/fwdpy11/pull/279>`_
* :class:`fwdpy11.TreeIterator` now allows iteration over sites and mutations in the current tree via :func:`fwdpy11.TreeIterator.sites` and :func:`fwdpy11.TreeIterator.mutations`. `PR 275 <https://github.com/molpopgen/fwdpy11/pull/275>`_
* Preorder traversal of nodes in the current tree is possible via :func:`fwdpy11.TreeIterator.nodes`.  Added :func:`fwdpy11.TreeIterator.samples` and :func:`fwdpy11.TreeIterator.samples_below`. `PR 272 <https://github.com/molpopgen/fwdpy11/pull/272>`_

0.5.0
+++++++++++

This is an intermediate release as we are still working towards supporting more general demographic models.

Major changes include:

* Updating the fwdpp back-end to the pre-release code for fwdpp 0.8.0.  Almost none of these changes are "user facing".
* Add :class:`fwdpy11.SiteTable`, :class:`fwdpy11.Site` and new fields to :class:`fwdpy11.MutationRecord`. `PR 258 <https://github.com/molpopgen/fwdpy11/pull/258>`_  These changes affect the API for some function calls. See :ref:`upgrade_path` for details.

Even though this release changes some of the tree sequence data structures, we are still able to read in files generated by version 0.4.5! (This is actually unit tested.)

Minor changes include:

* Add :func:`fwdpy11.gsl_version`. `PR 256 <https://github.com/molpopgen/fwdpy11/pull/256>`_
* :attr:`fwdpy11.Mutation.g` is converted to the mutation's age when dumping table collections to tskit's format. `PR 257 <https://github.com/molpopgen/fwdpy11/pull/257>`_
* New exception types from fwdpp registered as Python exceptions. `PR 260 <https://github.com/molpopgen/fwdpy11/pull/260>`_
* Several updates to documentation and to continuous integration testing.

0.4.5
+++++++++++

* :class:`fwdpy11.DataMatrixIterator` now correctly handles nested window coordinates. `PR 244 <https://github.com/molpopgen/fwdpy11/pull/244>`_.


0.4.4
+++++++++++

* Add :class:`fwdpy11.DataMatrixIterator`. `PR 243 <https://github.com/molpopgen/fwdpy11/pull/243>`_.
* Reduce time needed to execute unit tests of tree sequence functions.

0.4.3
++++++++++++++++++++++++++++++++++

* Minor fixes to packaging of source distrubition.
* Add a YCM config file to source repo
* Allow mutation and recombination regions to be empty. `PR 239 <https://github.com/molpopgen/fwdpy11/pull/239>`_.

0.4.2
++++++++++++++++++++++++++++++++++

Minor release:

* :class:`fwdpy11.VariantIterator`  may now skip neutral or selected sites during iteration. The behavior is specified
  by parameters passed to the class upon construction.
* Documentation updates

0.4.1
++++++++++++++++++++++++++++++++++

Minor release:

* Added position ranges to tree traversal.  `PR 232 <https://github.com/molpopgen/fwdpy11/pull/232>`_.
* Changed default type for range arguments for VariantIterator and data matrix generation. `PR 233 <https://github.com/molpopgen/fwdpy11/pull/233>`_.
* Skipping fixations is now optional in :func:`fwdpy11.data_matrix_from_tables`.
* The C++ back-end for population classes was changed to avoid deleting move constructors. `PR 231 <https://github.com/molpopgen/fwdpy11/pull/231>`_.
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

* :class:`fwdpy11.MlocusPop` is *tentatively* deprecated.  The new features described in :ref:`generalized_maps` make
  this class obsolete, but we will await a final verdict pending more testing.

Bug fixes
-------------------------------------------------------------

* A bug in handling fixations during simulations with tree sequence recording is fixed. This bug is 
  GitHub `Issue 200 <https://github.com/molpopgen/fwdpy11/issues/200>`_ and the fix is
  `PR 201 <https://github.com/molpopgen/fwdpy11/pull/201>`_.
* Updates to the fwdpp submodule fix a bug in :func:`fwdpy11.ts.infinite_sites`.  Previously, if the genome size 
  was not 1.0, then the number of mutations would be off by a factor of the genome size divided by 1.0.  The error was
  due to a bug upstream in fwdpp.
* A bug in how diploid metadata were updated by genetic value types has been fixed.  It is unlikely that this bug
  affected anyone unless they had written custom genetic value calculations where the offspring's genetic value 
  depended on the parental metadata. `PR 173 <https://github.com/molpopgen/fwdpy11/pull/173>`_. 

Support for multivariate mutational effects
-------------------------------------------------------------

`PR 164 <https://github.com/molpopgen/fwdpy11/pull/164>`_ introduced support for multidimensional mutational effects.
This pull request introduced several changes: 

The following new types are added:

* :class:`fwdpy11.MultivariateGaussianEffects`, which is a new "region" type
* :class:`fwdpy11.genetic_values.SlocusPopMultivariateGeneticValueWithMapping`, which is a new ABC for multivariate genetic values
* :class:`fwdpy11.genetic_values.MultivariateGeneticValueToFitnessMap`, which is a new ABC mapping multivariate trait values down to a (single) fitness value.
* :class:`fwdpy11.genetic_values.MultivariateGSS`, which is GSS based on the Euclidean distance from multiple optima
* :class:`fwdpy11.genetic_values.MultivariateGSSmo`, which is the multi-dimensional analog to the existing GSSmo
* :class:`fwdpy11.genetic_values.SlocusMultivariateEffectsStrictAdditive`, which is a new genetic value class for pleiotropic traits.

`PR 175 <https://github.com/molpopgen/fwdpy11/pull/175>`_ adds tracking of genetic values during simulation as numpy
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

* :ref:`pickling_pops`
* :ref:`generalized_maps`

Detailed changes:

* Add new function to pickle populations while using less memory. `PR 195 <https://github.com/molpopgen/fwdpy11/pull/195>`_,
  `PR 201 <https://github.com/molpopgen/fwdpy11/pull/201>`_
* Improved performance of simulations tracking lots of ancient samples. `PR 194 <https://github.com/molpopgen/fwdpy11/pull/194>`_
* Generalized genetic maps for single-locus simulations.  You can now do much of the "multi-locus" stuff with
  `SlocusPop` now. `PR 189 <https://github.com/molpopgen/fwdpy11/pull/189>`_
* Tree sequence recording now possible for mulit-locus simulations. `PR 185 <https://github.com/molpopgen/fwdpy11/pull/185>`_
* :func:`fwdpy11.ts.count_mutations` added. `PR 183 <https://github.com/molpopgen/fwdpy11/pull/183>`_, `PR 196 <https://github.com/molpopgen/fwdpy11/pull/196>`_, `PR 199 <https://github.com/molpopgen/fwdpy11/pull/199>`_
* Position and key properties added to :class:`fwdpy11.ts.VariantIterator`. `PR 180 <https://github.com/molpopgen/fwdpy11/pull/180>`_
  `PR 181 <https://github.com/molpopgen/fwdpy11/pull/181>`_
* :class:`fwdpy11.ts.TreeIterator` is added, which provides much faster tree traversal. `PR 176 <https://github.com/molpopgen/fwdpy11/pull/176>`_,
  `PR 177 <https://github.com/molpopgen/fwdpy11/pull/177>`_
* :func:`fwdpy11.ts.simplify` no longer retains ancient samples present in the input by default. To do so, explicitly
  label any ancient samples to retain as part of the the samples list passed to the function.
  `PR 169 <https://github.com/molpopgen/fwdpy11/pull/169>`_
* The types :class:`fwdpy11.Region` and :class:`fwdpy11.Sregion` have be re-implemented as C++-based classes, replacing 
  the previous pure Python classes.  `PR 163 <https://github.com/molpopgen/fwdpy11/pull/163>`_,
  `PR 174 <https://github.com/molpopgen/fwdpy11/pull/163>`_
* :attr:`fwdpy11.model_params.ModelParams.nregions` now defaults to an empty list, which simplifies setup for simulations
  with tree sequences. `commit <https://github.com/molpopgen/fwdpy11/commit/b557c4162cbfdfba6c9126ebec14c7f3f43884eb>`_. 
* When simulating with tree sequences, it is no longer an error to attempt to record ancient samples from the last
  generation of a simulation. `PR 162 <https://github.com/molpopgen/fwdpy11/pull/162>`_

Changes to the C++ back-end include:

* The genetic value types now store a vector of genetic values.  The idea is to generalize the type to handle both uni-
  and multi- variate genetic values. `PR 172 <https://github.com/molpopgen/fwdpy11/pull/172>`_

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
  modified. `PR #93 <https://github.com/molpopgen/fwdpy11/pull/93>`_
* Fix a bug in updating a population's mutation lookup table. This bug was upstream in fwdpp (`fwdpp issue 130 <https://github.com/molpopgen/fwdpp/issues/130>`_).  While definitely a bug, I could never find a case where simulation outputs were adversely affected.  In other words, simulation output remained the same after the fix, due to the rarity of the bug. `PR #98 <https://github.com/molpopgen/fwdpy11/pull/98>`_


API changes/new features:
----------------------------------------------------

* Added support for tree sequence recording.  `PR #142 <https://github.com/molpopgen/fwdpy11/pull/142>`_
* Populations may now be dumped/loaded to/from files. See :func:`fwdpy11.SlocusPop.dump_to_file` and
  :func:`fwdpy11.SlocusPop.load_from_file`.  Analagous functions exist for MlocusPop. `PR #148 <https://github.com/molpopgen/fwdpy11/pull/148>`_
* :func:`fwdpy11.SlocusPop.sample` and :func:`fwdpy11.MlocusPop.sample` now return a :class:`fwdpy11.sampling.DataMatrix`.
  `PR #118 <https://github.com/molpopgen/fwdpy11/pull/118>`_
* :class:`fwdpy11.sampling.DataMatrix` is refactored to match updates to fwdpp.  `PR #139 <https://github.com/molpopgen/fwdpy11/pull/139>`_
* :func:`fwdpy11.sampling.matrix_to_sample` now return a tuple with the neutral and selected data, respectively, as the
  two elements.  `PR #128 <https://github.com/molpopgen/fwdpy11/pull/128>`_
* Diploids have been refactored into two separate classes, :class:`fwdpy11.DiploidGenotype` and
  :class:`fwdpy11.DiploidMetadata`.  Both classes are valid NumPy dtypes.  See :ref:`processingpopsNP`. `PR #108 <https://github.com/molpopgen/fwdpy11/pull/108>`_
* :class:`fwdpy11.model_params.ModelParams` is massively simpilfied. There is now only one class! See :ref:`model_params`. `PR #108 <https://github.com/molpopgen/fwdpy11/pull/108>`_
* The design of objects related to calculating genetic values is vastly simplified.  See :ref:`genetic_values_types`. `PR #108 <https://github.com/molpopgen/fwdpy11/pull/108>`_
* Populations now contain functions to add mutations, replacing previous functions in fwdpy11.util.  `PR #94 <https://github.com/molpopgen/fwdpy11/pull/94>`_
* :class:`fwdpy11.MlocusPop` now requires that :attr:`fwdpy11.MlocusPop.locus_boundaries` be initialized upon
  construction. `PR #96 <https://github.com/molpopgen/fwdpy11/pull/96>`_
* The mutation position lookup table of a population is now a read-only property. See :ref:`mpos`. `PR #103 <https://github.com/molpopgen/fwdpy11/pull/103>`_
* The mutation position lookup table is now represented as a dict of lists. `PR #121 <https://github.com/molpopgen/fwdpy11/pull/121>`_
* A mutation or fixation can now be rapidy found by its "key".  See :func:`fwdpy11.Population.find_mutation_by_key`
  and :func:`fwdpy11.Population.find_fixation_by_key`.  `PR #106 <https://github.com/molpopgen/fwdpy11/pull/106>`_

Back-end changes
------------------------

* The build system now uses cmake.  `PR #151 <https://github.com/molpopgen/fwdpy11/pull/151>`_ `PR #152 <https://github.com/molpopgen/fwdpy11/pull/152>`_
* Most uses of C's assert macro are replaced with c++ exceptions.  `PR #141 <https://github.com/molpopgen/fwdpy11/pull/141>`_
* The C++ back-end of classes no longer contain any Python objects. `PR #114 <https://github.com/molpopgen/fwdpy11/pull/114>`_
* `PR #108 <https://github.com/molpopgen/fwdpy11/pull/108>`_ changes the back-end for representing diploids and for
  calculating genetic values.
* `PR #98 <https://github.com/molpopgen/fwdpy11/pull/98>`_ changes the definition of the populaton lookup table, using
  the same model as `fwdpp PR #132 <https://github.com/molpopgen/fwdpp/pull/132>`_
* Refactored class hierarchy for populations. `PR #85  <https://github.com/molpopgen/fwdpy11/pull/85>`_
* Updated to the fwdpp 0.6.x API and cleanup various messes that resulted. `PR #76 <https://github.com/molpopgen/fwdpy11/pull/76>`_ `PR #84 <https://github.com/molpopgen/fwdpy11/pull/84>`_ `PR #90 <https://github.com/molpopgen/fwdpy11/pull/90>`_ `PR #109 <https://github.com/molpopgen/fwdpy11/pull/109>`_ `PR #110 <https://github.com/molpopgen/fwdpy11/pull/110>`_
* The position of extinct variants is set to the max value of a C++ double. `PR #105 <https://github.com/molpopgen/fwdpy11/pull/105>`_
* An entirely new mutation type was introduced on the C++ side.  It is API compatible with the previous type (fwdpp's
  "popgenmut"), but has extra fields for extra flexibility. `PR #77 <https://github.com/molpopgen/fwdpy11/pull/77>`_ `PR #88 <https://github.com/molpopgen/fwdpy11/pull/88>`_
* Replaced `std::bind` with lambda closures for callbacks. `PR #80 <https://github.com/molpopgen/fwdpy11/pull/80>`_
* Fast exposure to raw C++ buffers improved for population objects. `PR #89 <https://github.com/molpopgen/fwdpy11/pull/89>`_
* Refactored long unit tests. `PR #91 <https://github.com/molpopgen/fwdpy11/pull/91>`_
* The GSL error handler is now turned off when fwdpy11 is imported and replaced with a custom handler to propagate GSL errors to C++ exceptions. `PR #140 <https://github.com/molpopgen/fwdpy11/pull/140>`_
* Population mutation position lookup table changed to an unordered multimap. `PR #102 <https://github.com/molpopgen/fwdpy11/pull/102>`_
* When a mutation is fixed or lost, its position is now set to the max value of a C++ double.  This change gets rid of
  some UI oddities when tracking mutations over time. `PR #106 <https://github.com/molpopgen/fwdpy11/pull/106>`_ and
  this `commit <https://github.com/molpopgen/fwdpy11/commit/96e8b6e7ca4b257cb8ae5e704f6a36a4b5bfa7bc>`_.

Version 0.1.4
++++++++++++++++++++++++++

Bug fixes:
--------------------------

* A bug affecting retrieval of multi-locus diploid key data as a buffer for numpy arrays is now fixed. `PR #72 <https://github.com/molpopgen/fwdpy11/pull/72>`_
* :attr:`fwdpy11.SingleLocusDiploid.label` is now pickled. `PR #34 <https://github.com/molpopgen/fwdpy11/pull/34>`_
    
API changes/new features:
----------------------------------------------------

* Population objects have new member functions ``sample`` and ``sample_ind``.  These replace
  :func:`fwdpy11.sampling.sample_separate`, which is now deprecated.  For example, see
  :func:`~fwdpy11.SlocusPop.sample` for more info. (The
  same member functions exist for *all* population objects.) `PR #62 <https://github.com/molpopgen/fwdpy11/pull/62>`_
* Improved support for pickling lower-level types. See the unit test file `tests/test_pickling.py` for examples of directly pickling things like mutations and containers of mutations.  `PR #55 <https://github.com/molpopgen/fwdpy11/pull/55>`_
* `__main__.py` added.  The main use is to help writing python modules based on fwdpy11. See :ref:`developers` for details. `PR #54 <https://github.com/molpopgen/fwdpy11/pull/54>`_
* Attributes `popdata` and `popdata_user` added to all population objects. `PR #52 <https://github.com/molpopgen/fwdpy11/pull/52>`_
* :attr:`fwdpy11.SingleLocusDiploid.parental_data` added as read-only field. `PR #51 <https://github.com/molpopgen/fwdpy11/pull/51>`_
* :attr:`fwdpy11.MlocusPop.locus_boundaries` is now writeable.
* :attr:`fwdpy11.sampling.DataMatrix.neutral` and :attr:`fwdpy11.sampling.DataMatrix.selected` are now writeable
  buffers. :attr:`fwdpy11.sampling.DataMatrix.ndim_neutral` and :attr:`fwdpy11.sampling.DataMatrix.ndim_selected` have
  been changed from functions to read-only properties. `PR #45 <https://github.com/molpopgen/fwdpy11/pull/45>`_
* The 'label' field of :class:`fwdpy11.Region` (and :class:`fwdpy11.Sregion`) now populate the label
  field of a mutation. `PR #32 <https://github.com/molpopgen/fwdpy11/pull/32>`_ See tests/test_mutation_labels.py for an example.
* Population objects may now be constructed programatically. See :ref:`popobjects`.   `PR #36 <https://github.com/molpopgen/fwdpy11/pull/36>`_ 

Back-end changes
------------------------

* The numpy dtype for :class:`fwdpy11.Mutation` has been refactored so that it generates tuples useable to construct object instances. This PR also removes some helper functions in favor of C++11 uniform initialization for these dtypes. `PR #72 <https://github.com/molpopgen/fwdpy11/pull/72>`_
* The documentation building process is greatly streamlined.  `PR #60 <https://github.com/molpopgen/fwdpy11/pull/60>`_
* Object namespaces have been refactored.  The big effect is to streamline the manual. `PR #59 <https://github.com/molpopgen/fwdpy11/pull/59>`_
* Travis CI now tests several Python versions using GCC 6 on Linux. `PR #44 <https://github.com/molpopgen/fwdpy11/pull/44>`_
* :func:`fwdpy11.wright_fisher_qtrait.evolve` has been updated to allow "standard popgen" models of multi-locus
  evolution. This change is a stepping stone to a future global simplification of the API. `PR #42 <https://github.com/molpopgen/fwdpy11/pull/42>`_
* The :class:`fwdpy11.Sregion` now store their callback data differently.  The result is a type that can be
  pickled in Python 3.6. `PR #39 <https://github.com/molpopgen/fwdpy11/pull/39>`_ 
* Travis builds are now Linux only and test many Python/GCC combos. `PR #38 <https://github.com/molpopgen/fwdpy11/pull/38>`_
* Update to fwdpp_ 0.5.7  `PR #35 <https://github.com/molpopgen/fwdpy11/pull/35>`_
* The method to keep fixations sorted has been updated so that the sorting is by position and fixation time. `PR #33 <https://github.com/molpopgen/fwdpy11/pull/33>`_
* The doctests are now run on Travis. `PR #30 <https://github.com/molpopgen/fwdpy11/pull/30>`_
* Removed all uses of placement new in favor of pybind11::pickle. `PR #26 <https://github.com/molpopgen/fwdpy11/pull/26>`_.
* fwdpy11 are now based on the @property/@foo.setter idiom for safety and code reuse.  `PR #21 <https://github.com/molpopgen/fwdpy11/pull/21>`_

Version 0.1.3.post1
++++++++++++++++++++++++++

* Fixed GitHub issues #23 and #25 via `PR #24 <https://github.com/molpopgen/fwdpy11/pull/24>`_.

Version 0.1.3
++++++++++++++++++++++++++

Bug fixes:
------------------------

* Issue #2 on GitHub fixed. [`commit <https://github.com/molpopgen/fwdpy11/commit/562a4d31947d9a7aae31f092ed8c014e94dc56db>`_]

API changes/new features:
------------------------------------------------

* :class:`fwdpy11.Sregion` may now model distrubitions of effect sizes on scales other than the effect size itself.  A scaling parameter allows the DFE to be functions of N, 2N, 4N, etc. [`PR #16 <https://github.com/molpopgen/fwdpy11/pull/16>`_]
  * Github issues 7, 8, and 9 resolved. All are relatively minor usability tweaks.
* :func:`fwdpy11.util.change_effect_size` added, allowing the "s" and "h" fields of :class:`fwdpy11.Mutation` to be changed. [`commit <https://github.com/molpopgen/fwdpy11/commit/ba4841e9407b3d98031801d7eea92b2661871eb2>`_].
* The attributes of :class:`fwdpy11.Mutation` are now read-only, addressing Issue #5 on GitHub. [`commit <https://github.com/molpopgen/fwdpy11/commit/f376d40788f3d59baa01d1d56b0aa99706560011>`_]
* Trait-to-fitness mapping functions for quantitative trait simulations now take the entire population, rather than just the generation.  This allows us to model things like truncation selection, etc. [`commit <https://github.com/molpopgen/fwdpy11/commit/fa37cb8f1763bc7f0e64c8620b6bc1ca350fddb9>`_]

Back-end changes
------------------------

* Code base updadted to work with pybind11_ 2.2.0. [`PR #19 <https://github.com/molpopgen/fwdpy11/pull/19>`_] 
* :mod:`fwdpy11.model_params` has been refactored, addressing issue #4 on GitHub.  The new code base is more idiomatic w.r.to Python's OO methods.`[`commit <https://github.com/molpopgen/fwdpy11/commit/1b811c33ab394ae4c64a3c8894984f320b870f22>`_]
* Many of the C++-based types can now be pickled, making model parameter objects easier to serialize.  Most of the
  changes are in [`this commit <https://github.com/molpopgen/fwdpy11/commit/d0a3602e71a866f7ff9d355d62953ea00c663c5a>`_].  This mostly addresses Issue #3 on GitHub.
* Added magic numbers to keep track of compatibility changes to serialization formats.
* __str__ changed to __repr__ for region types [`commit <https://github.com/molpopgen/fwdpy11/commit/2df859dd74d3de79d941a1cc21b8712a52bcf9ba>`_]
* fwdpy11.model_params now uses try/except rather than isinstance to check that rates are float-like types.[`commit <https://github.com/molpopgen/fwdpy11/commit/37112a60cd8fc74133945e522a47183314bf4085>`_]

Version 0.1.2
++++++++++++++++++++++++++

Bug fixes:
---------------------
* Fixed bug in setting the number of loci after deserializing a multi-locus population object. [`commit
  <https://github.com/molpopgen/fwdpy11/commit/4e4a547c5b4d30692b62bb4b4a5c22a4cd21d0fa>`_]

API and back-end changes:
------------------------------------------
* The C++ data structures are connected to NumPy via Python buffer protocol.  See :ref:`processingpopsNP`.  [`commit
  <https://github.com/molpopgen/fwdpy11/commit/48e3925a867c4ec55e1e5bb05457396fb456bc47>`_]
* :func:`fwdpy11.sampling.separate_samples_by_loci` changed to take a list of positions as first argument, and not a population object.

Version 0.1.1
++++++++++++++++++++++++++

Bug fixes:
---------------------
* Fixed bug in :func:`fwdpy11.sampling.DataMatrix.selected` that returned wrong data in best case scenario and could
  have caused crash in worst case. [`commit
  <https://github.com/molpopgen/fwdpy11/commit/e715fb74472555aa64e1d894563ec218ebba1a97>`_].
* Fix bug recording fixation times.  If a population was evolved multiple times, fixation times from the later rounds of
  evolution were incorrect. 
  [`commit <https://github.com/molpopgen/fwdpy11/commit/9db14d8b3db1c744045e20bfc00ce37e7fb28dfb>`_]
* Fix issue #1, related to fixations in quantitative trait sims. [`commit <https://github.com/molpopgen/fwdpy11/commit/6a27386498f056f0c4cc1fc6b8ea12f2b807636c>`_]
* The "label" field of a diploid is now initialized upon constructing a population.

API and back-end changes:
------------------------------------------
* Added :func:`fwdpy11.sampling.matrix_to_sample` and :func:`fwdpy11.sampling.separate_samples_by_loci`. [`commit <https://github.com/molpopgen/fwdpy11/commit/639c8de999679140fad6a976ff6c1996b25444aa>`_]
* Custom stateless fitness/genetic value calculations may now be implemented with a minimal amount of C++ code. See
  :ref:`customgvalues`. [`commit
  <https://github.com/molpopgen/fwdpy11/commit/a75166d9ff5471c2d18d66892f9fa01ebec5a667>`_]
* Custom fitness/genetic value calculations now allowed in pure Python, but they are quite slow (for now). See 
  :ref:`customgvalues`. [`commit <https://github.com/molpopgen/fwdpy11/commit/5549286046ead1181cba684464b3bcb19918321e>`_]
* Stateful trait value models enabled for qtrait sims. [`commit <https://github.com/molpopgen/fwdpy11/commit/161dfcef63f3abf28ad56df33b84a92d87d7750f>`_]
* Refactor evolution functions so that stateful fitness models behave as expected.  Enable compiling in a debug mode.
  Fix bug in operator== for diploid type. [`commit <https://github.com/molpopgen/fwdpy11/commit/a726c0535a5176aab1df5211fee7bf0aeba5054b>`_]
* fwdpy11.util added, providing :func:`fwdpy11.util.add_mutation`. [`commit <https://github.com/molpopgen/fwdpy11/commit/17b92dbe61ee85e2e60211e7dc0ed507a70dbd64>`_]
* Simulations now parameterized using classes in fwdpy11.model_params. [`commit <https://github.com/molpopgen/fwdpy11/commit/18e261c8596bf63d2d4e1ef228effb87397b793e>`_] and [`commit <https://github.com/molpopgen/fwdpy11/commit/eda7390adb9a98a5d96e6557ba1003488ebac511>`_]
* Added multi-locus simulation of quantitative traits. [`commit <https://github.com/molpopgen/fwdpy11/commit/fcad8de9d37bcef5a71ba6d26b4e40e1b67b1993>`_]
* Refactoring of type names. [`commit <https://github.com/molpopgen/fwdpy11/commit/632477c7b7592d956149a0cf44e4d26f2a67797e>`_]
* Refactoring internals of single-region fitness/trait value types. [`commit <https://github.com/molpopgen/fwdpy11/commit/d55d63631d02fdb2193940475dbcffaa201cf882>`_]
* Allow selected mutations to be retained in fwdpy11.wright_fisher.evolve_regions_sampler_fitness. [`commit <https://github.com/molpopgen/fwdpy11/commit/dcc1f2f6555eeada669efef8317f446e3cd0e46a>`_]

**Note:** the refactoring of type names will break scripts based on earlier versions.  Sorry, but things are rapidly changing here.  Please note that you can reassign class and function names in Python, allowing quick hacks to preserve compatibility:

.. code-block:: python

    import fwdpy11
    Spop = fwdpy11.SlocusPop

Alternately:

.. code-block:: python
    
    from fwdpy11 import SlocusPop as Spop

.. _pybind11: https://github.com/pybind/pybind11

