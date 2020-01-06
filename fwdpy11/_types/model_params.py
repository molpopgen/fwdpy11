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


class ModelParams(object):
    """
    Class to hold simulation parameters

    .. versionadded:: 0.1.1

    .. versionchanged:: 0.2.0
        Changed this from a horrible class hierarchy
        into a much simpler, single class.

    .. versionchanged:: 0.6.0

        Updated to support :class:`fwdpy11.DiscreteDemography`

    """

    def __init__(self, **kwargs):
        self.__nregions = []
        self.__sregions = None
        self.__recregions = None
        self.__demography = None
        self.__prune_selected = True
        self.__rates = None
        self.__gvalue = None
        self.__pself = 0.0
        self.__simlen = None
        self.__popsizes = None  # FIXME: this is a hack in 0.6.0
        for key, value in kwargs.items():
            if key in dir(self) and key[:1] != "_":
                setattr(self, key, value)
            elif key not in dir(self):
                raise ValueError(key, " not a valid parameter for this class.")

    def _set_demography(self, value):
        from fwdpy11 import DiscreteDemography
        if isinstance(value, DiscreteDemography):
            return value

        # Assume value is a numpy array of pop
        # sizes over time
        self.simlen = len(value)
        self.__popsizes = value  # FIXME: this is a hack in 0.6.0
        return DiscreteDemography(value)

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
        self.__demography = self._set_demography(value)

    @property
    def simlen(self):
        """
        Get or set the length of the simulation.

        .. versionadded:: 0.6.0
        """
        return self.__simlen

    @simlen.setter
    def simlen(self, value):
        self.__simlen = value

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
    def pself(self):
        """
        Get/set selfing probability/probabilities.
        """
        return self.__pself

    @pself.setter
    def pself(self, pself_):
        self.__pself = pself_

    @property
    def popsizes(self):
        """
        This is a hack in 0.6.0 so that
        simulations w/o tree sequences use the
        same API as previous versions.
        """
        # FIXME: hack introduced in 0.6.0
        return self.__popsizes

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
        if self.rates[1] > 0 and len(self.sregions) == 0:
            raise ValueError(
                "mutation rate to selected variants is > 0 but no Sregions are defined")
        if self.rates[0] > 0 and len(self.nregions) == 0:
            raise ValueError(
                "mutation rate to neutral variants is > 0 but no Regions are defined")
        if self.simlen is None or self.simlen < 0:
            raise ValueError("simlen cannot be none or < 0")
        # if self.rates[2] > 0 and len(self.recregions) == 0:
        #     raise ValueError("recombination rate is > 0 but no Regions are defined")
