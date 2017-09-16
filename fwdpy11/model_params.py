#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
#

"""
    Define classes for parameterizing simulations.

    The classes take either kwargs or an expanded dict as arguments.
    The names of the kwargs are the same as the names of class properties.
    The properties themselves are gettable/settable in the usual way.

    Each class has a validate function that attempts to do a
    full error-checking on the data stored in a class instance.

    The module fwdpy11.ezparams contains some functions for rapidly setting up
    dictionaries of parameters for common cases.

    .. todo::
        The classes should make deep copies of the input so that data
        from the calling environment are not changed, which can happen
        in some cases.
        There is an issue for C++ types holding
        things like Python objects.  Currently,
        tests/test_python_genetic_values.py
        will segfault if we copy.deepcopy the params objects.
"""

import warnings


def _validate_types(data, typename, strict):
    for i in data:
        if (strict is True and type(i) is not typename) or \
           isinstance(i, typename) is False:
            raise TypeError("invalid type: " + str(typename))


class ModelParams(object):
    """
    Base class for simulation parameters

    .. versionadded:: 0.1.1
    """

    def __init__(self, **kwargs):
        self.__nregions = None
        self.__sregions = None
        self.__recregions = None
        self.__demography = None
        self.__prune_selected = True
        for key, value in kwargs.items():
            if key in dir(self) and key[:1] != "_":
                setattr(self, key, value)
            elif key not in dir(self):
                raise ValueError(key, " not a valid parameter for this class.")

    @property
    def nregions(self):
        """
        Get or set the neutral regions.
        The implementation of the setter
        is handled by derived classes.
        Look there for more details.
        """
        return self.__nregions

    @property
    def sregions(self):
        """
        Get or set the selected regions.
        The implementation of the setter
        is handled by derived classes.
        Look there for more details.
        """
        return self.__sregions

    @property
    def recregions(self):
        """
        Get or set the recombination regions.
        The implementation of the setter
        is handled by derived classes.
        Look there for more details.
        """
        return self.__recregions

    @property
    def prune_selected(self):
        """
        Get or set whether or not selected fixations
        are to be removed during simulaton.

        When setting, a bool is required.
        """
        return self.__prune_selected

    @prune_selected.setter
    def prune_selected(self, value):
        self.__prune_selected = bool(value)

    @nregions.setter
    def nregions(self, nregions):
        self.__nregions = nregions

    @sregions.setter
    def sregions(self, sregions):
        self.__sregions = sregions

    @recregions.setter
    def recregions(self, recregions):
        self.__recregions = recregions

    @property
    def demography(self):
        """
        Get or set demographic history.
        The setters are fully implemented in
        derived classes.  Look there for
        type info, as the details are model-dependent.
        """
        return self.__demography

    @demography.setter
    def demography(self, value):
        self.__demography = value

    def validate(self):
        """
        Error check model params.

        :raises ValueError: Throws ValueError if validation fails.

        .. note:: Demography objects must be validated in derived classes.
        """
        if self.nregions is None:
            raise ValueError("neutral regions cannot be None")
        if self.sregions is None:
            raise ValueError("selected regions cannot be None")
        if self.recregions is None:
            raise ValueError("recombination regions cannot be None")
        if self.demography is None:
            raise ValueError("demography cannot be None")
        if self.prune_selected is None:
            raise ValueError("prune_selected cannot be None")


def _validate_single_deme_demography(value):
    import numpy as np
    if any(i < 0 for i in value):
        raise ValueError("all population sizes must be > 0")
    if (type(value) is np.ndarray) is False:
        raise ValueError("Type for population size " +
                         "history must be numpy.array")


def _validate_single_locus_rates(data):
    for key, value in data.items():
        if value is None:
            raise ValueError(key + " cannot be None")
        if value < 0.:
            raise ValueError("All mutation and recombination " +
                             "rates must be non-negative.")


class SlocusParams(ModelParams):
    """
    Simulation parameters for single-locus, single-deme simulations.

    .. note::
        The demography property for this class is a numpy
        array with dtype = np.uint32.

        When setting nregions or recregions, lists of
        instances of :class:`fwdpy11.regions.Region` are
        required.

        When setting sregions, lists of instances of
        :class:`fwdpy11.regions.Sregion` are required.

    .. note::
        If no gvalue is assigned, then
        :class:`fwdpy11.fitness.SlocusMult` with
        a scaling of 2.0 will be used.

    .. versionadded:: 0.1.1
    """

    def __init__(self, **kwargs):
        self.__mutrate_n = None
        self.__mutrate_s = None
        self.__recrate = None
        self.__gvalue = None
        self.__pself = 0.0
        gv_present = False
        if 'gvalue' in kwargs:
            gv_present = True
        super(SlocusParams, self).__init__(**kwargs)
        if gv_present is False:
            from fwdpy11.fitness import SlocusMult
            self.gvalue = SlocusMult(2.0)

    def __set_rates(self, input_rates):
        try:
            if len(input_rates) != 3:
                raise ValueError("len(input_rates) != 3")
            self.mutrate_n = input_rates[0]
            self.mutrate_s = input_rates[1]
            self.recrate = input_rates[2]
            _validate_single_locus_rates(self.rates)
        except:
            try:
                if isinstance(input_rates, dict):
                    for key, value in input_rates.items():
                        if key == "mutrate_n":
                            setattr(self, key, float(value))
                        elif key == "mutrate_s":
                            setattr(self, key, float(value))
                        elif key == "recrate":
                            setattr(self, key, float(value))
                        else:
                            raise ValueError("invalid key/value " +
                                             "pair for setting rates: ",
                                             key, ', ' + str(value))
                else:
                    raise
            except:
                raise
        _validate_single_locus_rates(self.rates)

    @property
    def mutrate_n(self):
        """
        Get or set the neutral mutation rate.
        When setting, a float >= 0.0 is required.
        """
        return self.__mutrate_n

    @mutrate_n.setter
    def mutrate_n(self, neutral_mutation_rate):
        try:
            float(neutral_mutation_rate)
        except:
            raise ValueError("neutral mutation rate must be a float")
        if neutral_mutation_rate < 0.0:
            raise ValueError("neutral_mutation_rate must be >= 0.0")
        self.__mutrate_n = neutral_mutation_rate

    @property
    def mutrate_s(self):
        """
        Get or set the selected mutation rate.
        When setting, a float >= 0.0 is required.
        """
        return self.__mutrate_s

    @mutrate_s.setter
    def mutrate_s(self, selected_mutation_rate):
        try:
            float(selected_mutation_rate)
        except:
            raise ValueError("selected mutation rate must be a float")
        if selected_mutation_rate < 0.0:
            raise ValueError("selected mutation rate must be >= 0.0")
        self.__mutrate_s = selected_mutation_rate

    @property
    def recrate(self):
        """
        Get or set the recombination rate.
        When setting, a float >= 0.0 is required.
        """
        return self.__recrate

    @recrate.setter
    def recrate(self, recombination_rate):
        try:
            float(recombination_rate)
        except:
            raise ValueError("recombination rate must be a float")
        if recombination_rate < 0.0:
            raise ValueError("recombination rate must be >= 0.0")
        self.__recrate = recombination_rate

    @property
    def pself(self):
        """
        Get or set selfing probability.
        When setting, a float in the interval
        [0,1] is required.
        """
        return self.__pself

    @pself.setter
    def pself(self, prob_selfing):
        try:
            float(prob_selfing)
        except:
            raise ValueError("selfing probability must be a float")
        if prob_selfing < 0.0 or prob_selfing > 1.0:
            raise ValueError(
                "invalid selfing probability: " + str(prob_selfing))
        self.__pself = float(prob_selfing)

    @property
    def rates(self):
        """
        Get or set all rates at once.
        When getting, a dict is returned.
        When setting, a list or tuple of the
        three rates (neutral, seleted, recombination)
        may be used, or a dict with the keys mutrate_n,
        mutrate_s, and recrate.
        """
        return {'mutrate_n': self.__mutrate_n,
                'mutrate_s': self.__mutrate_s,
                'recrate': self.__recrate}

    @rates.setter
    def rates(self, input_rates):
        self.__set_rates(input_rates)

    @property
    def gvalue(self):
        """
        Get or set genetic value calculator.
        When setting, an instance of
        :class:`fwdpy11.fitness.SlocusFitness`
        is required.
        """
        return self.__gvalue

    @gvalue.setter
    def gvalue(self, f):
        from fwdpy11.fitness import SlocusFitness
        if isinstance(f, SlocusFitness) is False:
            raise ValueError("invalid genetic value type: " + str(type(f)))
        self.__gvalue = f

    @ModelParams.demography.setter
    def demography(self, demog):
        _validate_single_deme_demography(demog)
        ModelParams.demography.fset(self, demog)

    @ModelParams.nregions.setter
    def nregions(self, neutral_regions):
        from fwdpy11 import Region
        _validate_types(neutral_regions, Region, True)
        ModelParams.nregions.fset(self, neutral_regions)

    @ModelParams.recregions.setter
    def recregions(self, recombination_regions):
        from fwdpy11 import Region
        _validate_types(recombination_regions, Region, True)
        ModelParams.recregions.fset(self, recombination_regions)

    @ModelParams.sregions.setter
    def sregions(self, selected_regions):
        from fwdpy11 import Sregion
        _validate_types(selected_regions, Sregion, False)
        ModelParams.sregions.fset(self, selected_regions)

    def validate(self):
        """
        Sanity-check parameters.
        """
        super(SlocusParams, self).validate()
        from fwdpy11.fitness import SlocusFitness

        _validate_single_locus_rates(self.rates)

        if callable(self.gvalue) is False:
            raise ValueError("genetic value function must be callable")

        # If individual rates are nonzero, then
        # corresponding lists of regions cannot be empty:
        checks = {'mutrate_n': self.nregions,
                  'mutrate_s': self.sregions,
                  'recrate': self.recregions}
        for key, value in self.rates.items():
            if value > 0.0 and len(checks[key]) == 0:
                raise ValueError(key + ' = ' + str(value) +
                                 ' is associated with empty list of regions.')
            # Issue warning as this may reflect a mistake:
            if value == 0. and len(checks[key]) > 0:
                warnings.warn(key + ' = ' + str(value) +
                              ' has nonzero region length ('
                              + str(checks[key]) +
                              ').  Is this intentional?')

        if isinstance(self.gvalue, SlocusFitness) is False:
            raise ValueError("invalid genetic value type: " +
                             type(self.__gvalue_data['gvalue']))

        _validate_single_deme_demography(self.demography)


class SlocusParamsQ(SlocusParams):
    """
    Single locus parameter object for quantitative trait simulations

    .. note::
        The demography property for this class is a numpy
        array with dtype = np.uint32.

        When setting nregions or recregions, lists of
        instances of :class:`fwdpy11.regions.Region` are
        required.

        When setting sregions, lists of instances of
        :class:`fwdpy11.regions.Sregion` are required.

    .. note::
        If no gvalue is assigned, then
        :class:`fwdpy11.trait_values.SlocusAdditiveTrait` with
        a scaling of 2.0 will be used.

        If no trait2w is assigned, then
        :class:`fwdpy11.wright_fisher_qtrait.GSS` with VS = 1.0
        and O = 0.0 will be used.

    .. versionadded:: 0.1.1
    """

    def __init__(self, **kwargs):
        self.__trait_to_fitness = None
        self.__noise = None
        gv_present = False
        if 'gvalue' in kwargs:
            gv_present = True
        trait2w_present = False
        if 'trait2w' in kwargs or 'trait_to_fitness' in kwargs:
            trait2w_present = True

        super(SlocusParamsQ, self).__init__(**kwargs)
        if gv_present is False:
            from fwdpy11.trait_values import SlocusAdditiveTrait
            self.gvalue = SlocusAdditiveTrait(2.0)
        if trait2w_present is False:
            from fwdpy11.wright_fisher_qtrait import GSS
            self.trait2w = GSS(VS=1.0, O=0.0)

        if self.prune_selected is True:
            warnings.warn("setting prune_selected to False", UserWarning)
            self.prune_selected = False

    @property
    def trait2w(self):
        """
        Get or set the trait value to fitness mapping function.
        When setting, a callable is required.
        """
        return self.__trait_to_fitness

    @property
    def noise(self):
        """
        Get or set the noise function.
        When setting, a callable is required.
        """
        return self.__noise

    @trait2w.setter
    def trait2w(self, trait_to_fitness):
        if callable(trait_to_fitness) is False:
            raise ValueError("trait_to_fitness must be a callable.")
        self.__trait_to_fitness = trait_to_fitness

    @property
    def trait_to_fitness(self):
        warnings.warn("\'trait_to_fitness\' property is deprecated. "
                      "Please use \'trait2w\'",
                      DeprecationWarning)
        return self.trait2w

    @trait_to_fitness.setter
    def trait_to_fitness(self, value):
        warnings.warn("\'trait_to_fitness\' property is deprecated. "
                      "Please use \'trait2w\'",
                      DeprecationWarning)
        self.trait2w = value

    @noise.setter
    def noise(self, noise):
        if callable(noise) is False:
            raise ValueError("noise must be a callable.")
        self.__noise = noise

    def validate(self):
        super(SlocusParamsQ, self).validate()

        if callable(self.trait2w) is False:
            raise ValueError("trait to fitness mapper must be a callable")

        if self.noise is not None and callable(self.noise) is False:
            raise ValueError("noise function must be callable")

        if self.prune_selected is True:
            raise ValueError("invalid value for prune_selected: " +
                             str(self.prune_selected))

        if self.__trait_to_fitness is None:
            raise ValueError("trait to fitness mapping " +
                             "function cannot be None")


def _validate_multilocus_rates(data):
    for key, value in data.items():
        if value is None:
            raise ValueError(key + " cannot be None")
        if isinstance(value, list) is False:
            raise ValueError(key + " must be a list of float")
        for i in value:
            try:
                float(i)
            except:
                raise ValueError(
                    "rates must be floating point values")
        if any(i < 0. for i in value):
            raise ValueError("all mutation and recombination rates " +
                             "must be non-negative.")
    # All rate lists must be equal in length
    lengths = set([len(value) for key, value in data.items()])
    if len(lengths) > 1:
        raise ValueError("all lists of mutation and recombination " +
                         "rates must be equal length.")


def _sanity_check_mutlilocus_rates_and_regions(rate_data,
                                               nregions,
                                               sregions,
                                               recregions):
    for i, j in zip(rate_data['mutrates_n'], nregions):
        if i > 0.0 and len(j) == 0:
            raise ValueError(
                "non-zero rate found associated with empty neutral locus")
    for i, j in zip(rate_data['mutrates_s'], sregions):
        if i > 0.0 and len(j) == 0:
            raise ValueError(
                "non-zero rate found associated with empty selected locus")
    for i, j in zip(rate_data['recrates'], recregions):
        if i > 0.0 and len(j) == 0:
            raise ValueError(
                "non-zero rate found associated with empty "
                "recombination locus")


class MlocusParams(ModelParams):
    """
    Model parameters for multi-locus simulation

    .. note::
        The demography property for this class is a numpy
        array with dtype = np.uint32.

        When setting nregions or recregions, lists of
        lists are required.  The inner lists may be empty
        or they may contain instances of
        :class:`fwdpy11.regions.Region`.

        Similarly, sregions inner lists should contain
        instances of :class:`fwdpy11.regions.Sregion`.
        They may also be empty.

    .. note::
        If aggregator is not assigned, then
        :class:`fwdpy11.multilocus.AggMultFitness` will be
        used.

        If gvalue is not assigned during construction, then
        :class:`fwdpy11.multilocus.MultiLocusGeneticValue`
        will be used and filled with instances of
        :class:`fwdpy11.fitness.SlocusMult` with
        scaling set to 2.0.  This assignment to gvalue is
        only possible if sregions is assigned during construction,
        so that we know how many loci there are.

    .. versionadded:: 0.1.1
    """

    def __init__(self, **kwargs):
        self.__mutrates_n = None
        self.__mutrates_s = None
        self.__recrates = None
        self.__interlocus = None
        self.__gvalue = None
        self.__agg = None
        self.__pself = 0.0
        super(MlocusParams, self).__init__(**kwargs)
        if self.aggregator is None:
            from fwdpy11.multilocus import AggMultFitness
            self.aggregator = AggMultFitness()
        if self.gvalue is None:
            if self.sregions is not None:
                from fwdpy11.multilocus import MultiLocusGeneticValue
                from fwdpy11.fitness import SlocusMult
                self.gvalue = MultiLocusGeneticValue([SlocusMult(2.0)] *
                                                     len(self.sregions))

    @property
    def gvalue(self):
        """
        Get or set genetic value calculation object.
        When setting, either a list of
        :class:`fwdpy11.fitness.SlocusFitness` or
        and instance of
        :class:`fwdpy11.multilocus.MultiLocusGeneticValue`
        is required.
        """
        return self.__gvalue

    @gvalue.setter
    def gvalue(self, f):
        from fwdpy11.multilocus import MultiLocusGeneticValue
        from fwdpy11.fitness import SlocusFitness
        if isinstance(f, list) is True:
            if any(isinstance(i, SlocusFitness) for i in f) is False:
                raise ValueError("all elements in list must " +
                                 "be single-locus genetic value objects.")
            self.__gvalue = MultiLocusGeneticValue(f)
        elif isinstance(f, MultiLocusGeneticValue) is True:
            self.__gvalue = f
        else:
            raise ValueError("invalid genetic value type: " + str(type(f)))

    @property
    def agg(self):
        """
        Get or set the aggregator function.
        When setting, a callable is required.

        .. note::
            This is deprecated, use 'aggregator'
            instead.

        .. deprecated:: 0.13.0
        """
        warnings.warn("\'agg\' property is deprecated. "
                      "Please use \'aggregator\'",
                      DeprecationWarning)
        return self.aggregator

    @agg.setter
    def agg(self, aggregator):
        warnings.warn("\'agg\' property is deprecated. "
                      "Please use \'aggregator\'",
                      DeprecationWarning)
        self.aggregator = aggregator

    @property
    def aggregator(self):
        """
        Get or set the aggregator function.
        When setting, a callable is required.
        """
        return self.__agg

    @aggregator.setter
    def aggregator(self, aggregator):
        if callable(aggregator) is False:
            raise ValueError("aggregator must be a callable")
        self.__agg = aggregator

    @property
    def mutrates_n(self):
        """
        Get or set the neutral mutation rates.
        When setting a list of floats >= 0.0
        is required.
        """
        return self.__mutrates_n

    @mutrates_n.setter
    def mutrates_n(self, neutral_mutation_rates):
        try:
            for i in neutral_mutation_rates:
                float(i)
        except:
            raise ValueError("neutral mutation rates must all be float")

        if any(i < 0.0 for i in neutral_mutation_rates):
            raise ValueError("neutral mutation rates must be >= 0.0")

        self.__mutrates_n = neutral_mutation_rates

    @property
    def mutrates_s(self):
        """
        Get or set the selected mutation rates.
        When setting a list of floats >= 0.0
        is required.
        """
        return self.__mutrates_s

    @mutrates_s.setter
    def mutrates_s(self, selected_mutation_rates):
        try:
            for i in selected_mutation_rates:
                float(i)
        except:
            raise ValueError("selected mutation rates must all be float")

        if any(i < 0.0 for i in selected_mutation_rates):
            raise ValueError("selected mutation rates must be >= 0.0")

        self.__mutrates_s = selected_mutation_rates

    @property
    def recrates(self):
        """
        Get or set the recombination rates.
        When setting a list of floats >= 0.0
        is required.
        """
        return self.__recrates

    @recrates.setter
    def recrates(self, recombination_rates):
        try:
            for i in recombination_rates:
                float(i)
        except:
            raise ValueError("recombination rates must all be float")

        if any(i < 0.0 for i in recombination_rates):
            raise ValueError("recombination rates must be >= 0.0")

        self.__recrates = recombination_rates

    @ModelParams.demography.setter
    def demography(self, demog):
        _validate_single_deme_demography(demog)
        ModelParams.demography.fset(self, demog)

    @ModelParams.nregions.setter
    def nregions(self, neutral_regions):
        from fwdpy11 import Region
        for i in neutral_regions:
            _validate_types(i, Region, True)
        ModelParams.nregions.fset(self, neutral_regions)

    @ModelParams.recregions.setter
    def recregions(self, recombination_regions):
        from fwdpy11 import Region
        for i in recombination_regions:
            _validate_types(i, Region, True)
        ModelParams.recregions.fset(self, recombination_regions)

    @ModelParams.sregions.setter
    def sregions(self, selected_regions):
        from fwdpy11 import Sregion
        for i in selected_regions:
            _validate_types(i, Sregion, False)
        ModelParams.sregions.fset(self, selected_regions)

    @property
    def pself(self):
        """
        Get or set selfing probability.
        When setting, a float in the interval
        [0,1] is required.
        """
        return self.__pself

    @pself.setter
    def pself(self, prob_selfing):
        self.__pself = float(prob_selfing)

    @property
    def rates(self):
        """
        Get or set all rates at once.
        When getting, a dict is returned.
        When setting, a list or tuple of the
        three rates (neutral, seleted, recombination)
        may be used, or a dict with the keys mutrates_n,
        mutrates_s, and recrates.
        """
        return {'mutrates_n': self.__mutrates_n,
                'mutrates_s': self.__mutrates_s,
                'recrates': self.__recrates}

    @rates.setter
    def rates(self, rates):
        self.__set_rates(rates)

    def __set_rates(self, rates):
        try:
            self.__mutrates_n = rates[0]
            self.__mutrates_s = rates[1]
            self.__recrates = rates[2]
            _validate_multilocus_rates(self.rates)
        except:
            try:
                if isinstance(rates, dict):
                    for key, value in rates.items():
                        if key in self.__expected_mutrec_kwargs:
                            self.__mutrec_data[key] = value
                        else:
                            raise ValueError("invalid key/value pair " +
                                             "for setting rates: ",
                                             key, ', ' + str(value))
                else:
                    raise
            except:
                raise
        _validate_multilocus_rates(self.__mutrec_data)

    @property
    def interlocus(self):
        """
        Get or set interlocus recombination rates.
        When setting, a list of instances of
        :class:`fwdpy11.multilocus.InterlocusRecombination`
        is required.
        """
        return self.__interlocus

    @interlocus.setter
    def interlocus(self, interlocus_rec):
        from fwdpy11.multilocus import InterlocusRecombination as IR
        if any(isinstance(i, IR) for i in interlocus_rec) is False:
            raise ValueError("interlocus recombination objects must be "
                             "instances of "
                             "fwdpy11.multilocus.InterlocusRecombination")
        self.__interlocus = interlocus_rec

    def validate(self):
        """
        Sanity-check the class data
        """
        super(MlocusParams, self).validate()

        _validate_multilocus_rates(self.rates)

        if callable(self.aggregator) is False:
            raise ValueError("aggregator must be callable")

        # Lengths of all regions must be the same
        rlens = set((len(self.nregions),
                     len(self.sregions), len(self.recregions)))
        if len(rlens) > 1:
            raise ValueError("length of all region " +
                             "containers must be equal: " +
                             ','.join([str(i) for i in rlens]))

        # Lengths of all lists of rates must be the same
        ratelens = set([len(value) for key, value in self.rates.items()])
        if len(ratelens) > 1:
            raise ValueError("length of per-locus " +
                             "rate containers must be equal: " +
                             ','.join([str(i) for i in rlens]))

        if rlens != ratelens:
            raise ValueError("region and rate containers of different lengths")

        _sanity_check_mutlilocus_rates_and_regions(self.rates,
                                                   self.nregions,
                                                   self.sregions,
                                                   self.recregions)

        # For k loci, there must be k-1 interlocus-recombination rates
        nloci = rlens.pop()
        if len(self.interlocus) + 1 != nloci:
            raise ValueError("invalid parameters: " +
                             str(nloci) + " loci, but only " +
                             str(len(self.interlocus)) +
                             " between-locus recombination functions.")

        if len(self.gvalue) != nloci:
            raise ValueError("Number of loci does not equal " +
                             "number of functions in genetic value object: " +
                             str(nloci) + " != " + str(len(self.aggregator)))


class MlocusParamsQ(MlocusParams):
    """
    Multi-locus model parameters for quantitative trait simulations

    .. note::
        The demography property for this class is a numpy
        array with dtype = np.uint32.

    .. note::
        If aggregator is not assigned, then
        :class:`fwdpy11.multilocus.AggMultTrait` will be
        used.

        If gvalue is not assigned during construction, then
        :class:`fwdpy11.multilocus.MultiLocusGeneticValue`
        will be used and filled with instances of
        :class:`fwdpy11.trait_values.SlocusAdditiveTrait` with
        scaling set to 2.0.  This assignment to gvalue is
        only possible if sregions is assigned during construction,
        so that we know how many loci there are.

        If trait2w is not set during construction, then
        :class:`fwdpy11.wright_fisher_qtrait.GSS` will be
        used with VS = 1.0 and O = 0.0.

    .. versionadded:: 0.1.1
    """

    def __init__(self, **kwargs):
        self.__trait2w = None
        self.__noise = None
        agg_present = False
        if 'agg' in kwargs or 'aggregator' in kwargs:
            agg_present = True
        gvalue_present = False
        if 'gvalue' in kwargs:
            gvalue_present = True
        trait2w_present = False
        if 'trait2w' in kwargs or 'trait_to_fitness' in kwargs:
            trait2w_present = True
        super(MlocusParamsQ, self).__init__(**kwargs)
        if agg_present is False:
            from fwdpy11.multilocus import AggMultTrait
            self.aggregator = AggMultTrait()
        if gvalue_present is False and self.sregions is not None:
            from fwdpy11.multilocus import MultiLocusGeneticValue
            from fwdpy11.trait_values import SlocusAdditiveTrait
            self.gvalue = MultiLocusGeneticValue([SlocusAdditiveTrait(2.0)] *
                                                 len(self.sregions))
        if trait2w_present is False:
            from fwdpy11.wright_fisher_qtrait import GSS
            self.trait2w = GSS(VS=1.0, O=0.0)

    @property
    def trait2w(self):
        """
        Get or set the trait value to fitness mapping function.
        When setting, a callable is required.
        """
        return self.__trait2w

    @trait2w.setter
    def trait2w(self, trait2w):
        if callable(trait2w) is False:
            raise ValueError("trait2w must be a callable.")
        self.__trait2w = trait2w

    @property
    def trait_to_fitness(self):
        """
        Get or set the trait value to fitness mapping function.
        When setting, a callable is required.

        .. deprecated:: 0.13.0
            Use trait2w instead.
        """
        warnings.warn("\'trait_to_fitness\' property is deprecated. "
                      "Please use \'trait2w\'",
                      DeprecationWarning)
        return self.trait2w

    @trait_to_fitness.setter
    def trait_to_fitness(self, value):
        warnings.warn("\'trait_to_fitness\' property is deprecated. "
                      "Please use \'trait2w\'",
                      DeprecationWarning)
        self.trait2w = value

    @property
    def noise(self):
        """
        Get or set the noise function.
        When setting, a callable is required.
        """
        return self.__noise

    @noise.setter
    def noise(self, noise):
        if callable(noise) is False:
            raise ValueError("noise must be a callable.")
        self.__noise = noise

    def validate(self):
        """
        Sanity-check class data
        """
        if callable(self.trait2w) is False:
            raise ValueError("trait2w must be a callable.")

        if self.noise is None:
            warnings.warn("noise parameter is None.  Defaults will be used.",
                          UserWarning)

        super(MlocusParamsQ, self).validate()
