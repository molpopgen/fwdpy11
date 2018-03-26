Changelog
====================================================================================

Major changes are listed below.  Each release likely contains fiddling with back-end code, updates to latest fwdpp
version, etc.

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
* Added :func:`fwdpy11.sampling.matrix_to_sample` and :func:`fwdpy11.sampling.separate_samples_by_loci`. [`commit <https://github.com/molpopgen/fwdpy11/commit/i639c8de999679140fad6a976ff6c1996b25444aa>`_]
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
