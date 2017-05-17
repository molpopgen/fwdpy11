Changelog
====================================================================================

Major changes are listed below.  Each release likely contains fiddling with back-end code, updates to latest fwdpp
version, etc.

Version 0.1.1
++++++++++++++++++++++++++

Bug fixes:
---------------------
* Fix bug recording fixation times.  If a population was evolved multiple times, fixation times from the later rounds of
  evolution were incorrect. 
  [`commit <https://github.com/molpopgen/fwdpy11/commit/9db14d8b3db1c744045e20bfc00ce37e7fb28dfb>`_]
* Fix issue #1, related to fixations in quantitative trait sims. [`commit <https://github.com/molpopgen/fwdpy11/commit/6a27386498f056f0c4cc1fc6b8ea12f2b807636c>`_]
* The "label" field of a diploid is now initialized upon constructing a population.

API and back-end changes:
------------------------------------------
* Custom fitness/genetic value calculations now allowed in pure Python, but they are quite slow (for now). See 
  :ref:`customgvalue`. [`commit <https://github.com/molpopgen/fwdpy11/commit/5549286046ead1181cba684464b3bcb19918321e>`]_
* Stateful trait value models enabled for qtrait sims. [`commit <https://github.com/molpopgen/fwdpy11/commit/161dfcef63f3abf28ad56df33b84a92d87d7750f>`_]
* Refactor evolution functions so that stateful fitness models behave as expected.  Enable compiling in a debug mode.
  Fix bug in operator== for diploid type. [`commit <https://github.com/molpopgen/fwdpy11/commit/a726c0535a5176aab1df5211fee7bf0aeba5054b>`_]
* fwdpy11.util added, providing :func:`fwdpy11.util.add_mutation`. [`commit <https://github.com/molpopgen/fwdpy11/commit/17b92dbe61ee85e2e60211e7dc0ed507a70dbd64>`_]
* Simulations now parameterized using classes in fwdpy11.model_params. [`commit <https://github.com/molpopgen/fwdpy11/commit/18e261c8596bf63d2d4e1ef228effb87397b793e>`_] and [`commit <https://github.com/molpopgen/fwdpy11/commit/eda7390adb9a98a5d96e6557ba1003488ebac511>`_]
* Added multi-locus simulation of quantitative traits. [`commit <https://github.com/molpopgen/fwdpy11/commit/fcad8de9d37bcef5a71ba6d26b4e40e1b67b1993>`_]
* Refactoring of type names. [`commit <https://github.com/molpopgen/fwdpy11/commit/632477c7b7592d956149a0cf44e4d26f2a67797e>`_]
* Refactoring internals of single-region fitnes/trait value typess. [`commit <https://github.com/molpopgen/fwdpy11/commit/d55d63631d02fdb2193940475dbcffaa201cf882>`_]
* Allow selected mutations to be retained in fwdpy11.wright_fisher.evolve_regions_sampler_fitness. [`commit <https://github.com/molpopgen/fwdpy11/commit/dcc1f2f6555eeada669efef8317f446e3cd0e46a>`_]

**Note:** the refactoring of type names will break scripts based on earlier versions.  Sorry, but things are rapidly changing here.  Please note that you can reassign class and function names in Python, allowing quick hacks to preserve compatibility:

.. code-block:: python

    import fwdpy11
    Spop = fwdpy11.SlocusPop

Alternately:

.. code-block:: python
    
    from fwdpy11 import SlocusPop as Spop
