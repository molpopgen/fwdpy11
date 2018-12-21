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


def _validate_types(data, typename, strict):
    for i in data:
        if (strict is True and type(i) is not typename) or \
           isinstance(i, typename) is False:
            raise TypeError("invalid type: " + str(typename))


class ModelParams(object):
    """
    Class to hold simulation parameters

    .. versionadded:: 0.1.1

    .. versionchanged:: 0.2.0
        Changed this from a horrible class hierarchy
        into a much simpler, single class.

    See :ref:`model_params` for detailed documentation.
    """

    def __init__(self, **kwargs):
        self.__nregions = []
        self.__sregions = None
        self.__recregions = None
        self.__interlocus_rec = None
        self.__demography = None
        self.__prune_selected = True
        self.__rates = None
        self.__gvalue = None
        self.__pself = 0.0
        for key, value in kwargs.items():
            if key in dir(self) and key[:1] != "_":
                setattr(self, key, value)
            elif key not in dir(self):
                raise ValueError(key, " not a valid parameter for this class.")

    @property
    def nregions(self):
        """
        Get or set the neutral regions.
        """
        return self.__nregions

    @nregions.setter
    def nregions(self, nregions_):
        self.__nregions = nregions_

    @property
    def sregions(self):
        """
        Get or set the selected regions.
        """
        return self.__sregions

    @sregions.setter
    def sregions(self, sregions_):
        self.__sregions = sregions_

    @property
    def recregions(self):
        """
        Get or set the recombination regions.
        """
        return self.__recregions

    @recregions.setter
    def recregions(self, recregions_):
        self.__recregions = recregions_

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

    @property
    def demography(self):
        """
        Get or set demographic history.
        """
        return self.__demography

    @demography.setter
    def demography(self, value):
        self.__demography = value

    @property
    def gvalue(self):
        """
        Get/set the type name of the genetic value
        calculator.
        """
        return self.__gvalue

    @gvalue.setter
    def gvalue(self, gvalue_):
        self.__gvalue = gvalue_

    @property
    def rates(self):
        """
        Get/set mutation and recombination rates.
        """
        return self.__rates

    @rates.setter
    def rates(self, rates_):
        if len(rates_) != 3:
            raise ValueError("length of rates must be 3")
        self.__rates = rates_

    @property
    def mutrate_n(self):
        """
        Return neutral mutation rate(s).
        """
        return self.__rates[0]

    @property
    def mutrate_s(self):
        """
        Return selected mutation rate(s).
        """
        return self.__rates[1]

    @property
    def recrate(self):
        """
        Return recombination rate(s).
        """
        return self.__rates[2]

    @property
    def mutrates_n(self):
        """
        Return neutral mutation rate(s).
        """
        return self.mutrate_n

    @property
    def mutrates_s(self):
        """
        Return selected mutation rate(s).
        """
        return self.mutrate_s

    @property
    def recrates(self):
        """
        Return recombination rate(s).
        """
        return self.recrate

    @property
    def interlocus_rec(self):
        """
        Get/set interlocus recombination functions.
        """
        return self.__interlocus_rec

    @interlocus_rec.setter
    def interlocus_rec(self, ir):
        self.__interlocus_rec = ir

    @property
    def pself(self):
        """
        Get/set selfing probability/probabilities.
        """
        return self.__pself

    @pself.setter
    def pself(self, pself_):
        self.__pself = pself_

    def validate(self):
        """
        Error check model params.

        :raises TypeError: Throws TypeError if validation fails.

        """
        if self.nregions is None:
            raise TypeError("neutral regions cannot be None")
        if self.sregions is None:
            raise TypeError("selected regions cannot be None")
        if self.recregions is None:
            raise TypeError("recombination regions cannot be None")
        if self.demography is None:
            raise TypeError("demography cannot be None")
        if self.prune_selected is None:
            raise TypeError("prune_selected cannot be None")
        if self.gvalue is None:
            raise TypeError("gvalue cannot be None")
        if self.rates is None:
            raise TypeError("rates cannot be None")


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
