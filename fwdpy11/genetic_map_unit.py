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

import attr
import numpy as np

import fwdpy11._fwdpy11

from .class_decorators import attr_class_to_from_dict

_common_attr_attribs = {"frozen": True, "auto_attribs": True, "repr_ns": "fwdpy11"}


@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class PoissonInterval(fwdpy11._fwdpy11._ll_PoissonInterval):
    """
    Generate poisson number of crossover breakpoints.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    mean: float

    def __attrs_post_init__(self):
        super(PoissonInterval, self).__init__(self.beg, self.end, self.mean)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(PoissonInterval, self).__init__(**d)


@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class PoissonPoint(fwdpy11._fwdpy11._ll_PoissonPoint):
    """
    Generate a recombination breakpoint at a fixed position if the
    number of crossover events is odd.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    position: float
    mean: float

    def __attrs_post_init__(self):
        super(PoissonPoint, self).__init__(self.position, self.mean)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(PoissonPoint, self).__init__(**d)


@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class BinomialInterval(fwdpy11._fwdpy11._ll_BinomialInterval):
    """
    Generate exactly one crossover with a given probability

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    .. versionadded:: 0.5.2

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    probability: float

    def __attrs_post_init__(self):
        super(BinomialInterval, self).__init__(self.beg, self.end, self.probability)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(BinomialInterval, self).__init__(**d)


@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class BinomialPoint(fwdpy11._fwdpy11._ll_BinomialPoint):
    """
    Generate a crossover breakpoint at a fixed position with a
    fixed probability.  This class represents genetic distance
    as centiMorgans/100.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

            Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    position: float
    probability: float

    def __attrs_post_init__(self):
        super(BinomialPoint, self).__init__(self.position, self.probability)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(BinomialPoint, self).__init__(**d)


@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class FixedCrossovers(fwdpy11._fwdpy11._ll_FixedCrossovers):
    """
    Generate a fixed number of crossover breakpoints.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.5.0

        Refactored back-end to be based on fwdpp types

    .. versionchanged:: 0.7.1

        Add __repr__

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    num_xovers: int

    def __attrs_post_init__(self):
        super(FixedCrossovers, self).__init__(self.beg, self.end, self.num_xovers)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        super(FixedCrossovers, self).__init__(**d)
