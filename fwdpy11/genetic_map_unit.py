#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

import typing

import attr
import numpy as np

import fwdpy11._fwdpy11

from .class_decorators import (attr_add_asblack, attr_class_pickle_with_super,
                               attr_class_to_from_dict)


def _is_integer_if_discrete(self, attribute, value):
    if self.discrete is True:
        attr.validators.instance_of(int)(self, attribute, value)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class PoissonInterval(fwdpy11._fwdpy11._ll_PoissonInterval):
    """
    Generate poisson number of crossover breakpoints.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: The beginning of the region
    :type beg: int or float
    :param end: The end of the region
    :type end: int or float
    :param mean: The mean number of breakpoints per meiosis
    :type mean: float
    :param discrete: If `False`, positions are continuous and uniform from `[beg, end)`.
     If `True`, positions take integer values uniformly from `[beg, end)`.
    :type discrete: bool

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class

    .. versionchanged:: 0.12.0

        Added `discrete` option to initializer.
    """

    beg: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    end: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    mean: float
    discrete: bool = attr.ib(kw_only=True, default=False)

    def __attrs_post_init__(self):
        super(PoissonInterval, self).__init__(
            beg=self.beg, end=self.end, mean=self.mean, discrete=self.discrete
        )


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class PoissonPoint(fwdpy11._fwdpy11._ll_PoissonPoint):
    """
    Generate a recombination breakpoint at a fixed position if the
    number of crossover events is odd.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param position: The position of the crossover
    :type position: int or float
    :param mean: The mean number of breakpoints per meiosis
    :type mean: float
    :param discrete: If `False`, positions are continuous and uniform from `[beg, end)`.
     If `True`, positions take integer values uniformly from `[beg, end)`.
    :type discrete: bool

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class

    .. versionchanged:: 0.12.0

        Added `discrete` option to initializer.
    """

    position: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    mean: float
    discrete: bool = attr.ib(kw_only=True, default=False)

    def __attrs_post_init__(self):
        super(PoissonPoint, self).__init__(
            position=self.position, mean=self.mean, discrete=self.discrete
        )

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(PoissonPoint, self).__init__(**d)


@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class BinomialInterval(fwdpy11._fwdpy11._ll_BinomialInterval):
    """
    Generate exactly one crossover with a given probability

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: The beginning of the region
    :type beg: int or float
    :param end: The end of the region
    :type end: int or float
    :param probability: The probability of a recombination (per meiosis).
    :type probability: float
    :param discrete: If `False`, positions are continuous and uniform from `[beg, end)`.
     If `True`, positions take integer values uniformly from `[beg, end)`.
    :type discrete: bool

    .. versionadded:: 0.5.2

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class

    .. versionchanged:: 0.12.0

        Added `discrete` option to initializer.
    """

    beg: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    end: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    probability: float
    discrete: bool = attr.ib(kw_only=True, default=False)

    def __attrs_post_init__(self):
        super(BinomialInterval, self).__init__(
            beg=self.beg,
            end=self.end,
            probability=self.probability,
            discrete=self.discrete,
        )


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class BinomialPoint(fwdpy11._fwdpy11._ll_BinomialPoint):
    """
    Generate a crossover breakpoint at a fixed position with a
    fixed probability.  This class represents genetic distance
    as centiMorgans/100.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param position: The beginning of the region
    :type position: int or float
    :param probability: The probability of a recombination (per meiosis).
    :type probability: float
    :param discrete: If `False`, positions are continuous and uniform from `[beg, end)`.
     If `True`, positions take integer values uniformly from `[beg, end)`.
    :type discrete: bool

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

            Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class

    .. versionchanged:: 0.12.0

        Added `discrete` option to initializer.
    """

    position: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    probability: typing.Union[int, float]
    discrete: bool = attr.ib(kw_only=True, default=False)

    def __attrs_post_init__(self):
        super(BinomialPoint, self).__init__(
            position=self.position, probability=self.probability, discrete=self.discrete
        )


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class FixedCrossovers(fwdpy11._fwdpy11._ll_FixedCrossovers):
    """
    Generate a fixed number of crossover breakpoints.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: The beginning of the region
    :type beg: int or float
    :param end: The end of the region
    :type end: int or float
    :param num_xovers: The number of breakpoints per meiosis
    :type num_xovers: float
    :param discrete: If `False`, positions are continuous and uniform from `[beg, end)`.
     If `True`, positions take integer values uniformly from `[beg, end)`.
    :type discrete: bool

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class

    .. versionchanged:: 0.12.0

        Added `discrete` option to initializer.
    """

    beg: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    end: typing.Union[int, float] = attr.ib(validator=_is_integer_if_discrete)
    num_xovers: int
    discrete: bool = attr.ib(kw_only=True, default=False)

    def __attrs_post_init__(self):
        super(FixedCrossovers, self).__init__(
            beg=self.beg,
            end=self.end,
            num_xovers=self.num_xovers,
            discrete=self.discrete,
        )
