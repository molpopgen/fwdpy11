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
import warnings

import attr
import numpy as np

import fwdpy11._fwdpy11

from .class_decorators import (attr_add_asblack, attr_class_pickle_with_super,
                               attr_class_to_from_dict,
                               attr_class_to_from_dict_no_recurse)


def _add_deprecated_properties(self):
    """
    We've been using some non-Pythonic property names for a while.
    This decorator adds them back
    """

    def b(self):
        warnings.warn(
            "b is deprecated in favor of beg and is scheduled for removal in 0.11.",
            FutureWarning,
        )
        return self.beg

    def e(self):
        warnings.warn(
            "e is deprecated in favor of end and is scheduled for removal in 0.11.",
            FutureWarning,
        )
        return self.end

    def w(self):
        warnings.warn(
            "w is deprecated in favor of weight and is scheduled for removal in 0.11.",
            FutureWarning,
        )
        return self.weight

    def c(self):
        warnings.warn(
            "c is deprecated in favor of coupled and is scheduled for removal in 0.11.",
            FutureWarning,
        )
        return self.coupled

    def l(self):
        warnings.warn(
            "l is deprecated in favor of label and is scheduled for removal in 0.11.",
            FutureWarning,
        )
        return self.label

    self.b = property(b)
    self.e = property(e)
    self.w = property(w)
    self.c = property(c)
    self.l = property(l)

    return self


_common_attr_attribs = {"frozen": True, "auto_attribs": True, "repr_ns": "fwdpy11"}


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class Region(fwdpy11._fwdpy11._ll_Region):
    """
    A genomic region, defined by half-open interval [beg, end)

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: the beginning of the region
    :type beg: float
    :param end: the end of the region
    :type end: float
    :param weight: the weight to assign
    :type weight: float
    :param coupled: if True, the weight is converted to (end-beg)*weight
    :type coupled: bool
    :param label: Not relevant to recombining regions.
                  Otherwise, this value will be used to fill
                  :attr:`fwdpy11.Mutation.label`.
    :type label: numpy.uint16

    When coupled is True, the "weight" may be interpreted
    as a "per base pair" (or per unit, generally speaking) term.

    .. versionchanged:: 0.3.0

        Refactored from a pure Python class to a C++/pybind11 class

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    weight: float
    coupled: bool = attr.ib(default=True)
    label: int = attr.ib(default=0)

    def __attrs_post_init__(self):
        super(Region, self).__init__(
            self.beg, self.end, self.weight, self.coupled, self.label
        )


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class ConstantS(fwdpy11._fwdpy11._ll_ConstantS):
    """
    Mutations with fixed effect sizes

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: the beginning of the region
    :type beg: float
    :param end: the end of the region
    :type end: float
    :param weight: the weight to assign
    :type weight: float
    :param s: the selection coefficient
    :type s: float
    :param h: the dominance
    :type h: float
    :param coupled: if True, the weight is converted to (end-beg)*weight
    :type coupled: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16
    :param scaling: The scaling of the DFE
    :type scaling: float

    When coupled is True, the "weight" may be interpreted
    as a "per base pair" (or per unit, generally speaking) term.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    weight: float
    s: float
    h: float = attr.ib(default=1.0)
    coupled: bool = attr.ib(default=True)
    label: int = attr.ib(default=0)
    scaling: float = attr.ib(default=1.0)

    def __attrs_post_init__(self):
        super(ConstantS, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.s,
            self.h,
            self.coupled,
            self.label,
            self.scaling,
        )


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class ExpS(fwdpy11._fwdpy11._ll_ExpS):
    """
    Exponential distribution of effect sizes

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: the beginning of the region
    :type beg: float
    :param end: the end of the region
    :type end: float
    :param weight: the weight to assign
    :type weight: float
    :param mean: the mean selection coefficient
    :type s: float
    :param h: the dominance
    :type h: float
    :param coupled: if True, the weight is converted to (end-beg)*weight
    :type coupled: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16
    :param scaling: The scaling of the DFE
    :type scaling: float

    When coupled is True, the "weight" may be interpreted
    as a "per base pair" (or per unit, generally speaking) term.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    import numpy as np

    beg: float
    end: float
    weight: float
    mean: float
    h: float = attr.ib(default=1.0)
    coupled: bool = attr.ib(default=True)
    label: int = attr.ib(default=0)
    scaling: float = attr.ib(default=1.0)

    def __attrs_post_init__(self):
        super(ExpS, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.mean,
            self.h,
            self.coupled,
            self.label,
            self.scaling,
        )


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class GammaS(fwdpy11._fwdpy11._ll_GammaS):
    """
    Gamma distribution of effect sizes

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: the beginning of the region
    :type beg: float
    :param end: the end of the region
    :type end: float
    :param weight: the weight to assign
    :type weight: float
    :param mean: the mean selection coefficient
    :type mean: float
    :param shape: the shape parameter of the distribution
    :type shape: float
    :param h: the dominance
    :type h: float
    :param coupled: if True, the weight is converted to (end-beg)*weight
    :type coupled: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16
    :param scaling: The scaling of the DFE
    :type scaling: float

    When coupled is True, the "weight" may be
    interpreted as a "per base pair" (or per unit, generally speaking) term.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    weight: float
    mean: float
    shape_parameter: float
    h: float = 1.0
    coupled: bool = True
    label: int = 0
    scaling: float = 1.0

    def __attrs_post_init__(self):
        super(GammaS, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.mean,
            self.shape_parameter,
            self.h,
            self.coupled,
            self.label,
            self.scaling,
        )


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class GaussianS(fwdpy11._fwdpy11._ll_GaussianS):
    """
    Gaussian distribution of effect sizes

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: the beginning of the region
    :type beg: float
    :param end: the end of the region
    :type end: float
    :param weight: the weight to assign
    :type weight: float
    :param sd: standard deviation of effect sizes
    :type sd: float
    :param h: the dominance
    :type h: float
    :param coupled: if True, the weight is converted to (end-beg)*weight
    :type coupled: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16
    :param scaling: The scaling of the DFE
    :type scaling: float

    When coupled is True, the "weight" may be
    interpreted as a "per base pair" (or per unit, generally speaking) term.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    weight: float
    sd: float
    h: float = 1.0
    coupled: bool = True
    label: int = 0
    scaling: float = 1.0

    def __attrs_post_init__(self):
        super(GaussianS, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.sd,
            self.h,
            self.coupled,
            self.label,
            self.scaling,
        )


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class LogNormalS(fwdpy11._fwdpy11._ll_LogNormalS):
    """
    Log-normal distribution of effect sizes.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: the beginning of the region
    :type beg: float
    :param end: the end of the region
    :type end: float
    :param weight: the weight to assign
    :type weight: float
    :param zeta: the zeta parameter
    :type zeta: float
    :param sigma: the sigma parameter
    :type sigma: float
    :param h: the dominance
    :type h: float
    :param coupled: if True, the weight is converted to(end-beg)*weight
    :type Constructor: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16
    :param scaling: The scaling of the DFE
    :type scaling: float

    When coupled is True, the "weight" may be
    interpreted as a "per base pair" (or per unit, generally speaking) term.

    .. versionadded:: 0.7.0

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    weight: float
    zeta: typing.Optional[float]
    sigma: typing.Optional[float]
    h: float = 1.0
    coupled: bool = True
    label: int = 0
    scaling: float = 1.0

    def __attrs_post_init__(self):
        super(LogNormalS, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.zeta,
            self.sigma,
            self.h,
            self.coupled,
            self.label,
            self.scaling,
        )

    @classmethod
    def mv(
        cls,
        beg: float,
        end: float,
        weight: float,
        h: float = 1.0,
        coupled: bool = True,
        label: int = 0,
        scaling: float = 1.0,
    ):
        """
        Factory method to create an instance compatible with
        :class:`fwdpy11.mvDES`. See :ref:`mvdes` for
        details.

        :param beg: the beginning of the region
        :type beg: float
        :param end: the end of the region
        :type end: float
        :param weight: the weight to assign
        :type weight: float
        :param h: the dominance
        :type h: float
        :param coupled: if True, the weight is converted to(end-beg)*weight
        :type coupled: bool
        :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
        :type label: numpy.uint16
        :param scaling: The scaling of the DFE
        :type scaling: float
        """

        return cls(
            beg=beg,
            end=end,
            weight=weight,
            zeta=None,
            sigma=None,
            h=h,
            coupled=coupled,
            label=label,
            scaling=scaling,
        )


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class UniformS(fwdpy11._fwdpy11._ll_UniformS):
    """
    Uniform distrubution of effect sizes

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: the beginning of the region
    :type beg: float
    :param end: the end of the region
    :type end: float
    :param weight: the weight to assign
    :type weight: float
    :param lo: lower bound on s
    :type lo: float
    :param hi: upper bound on s
    :type hi: float
    :param h: the dominance
    :type h: float
    :param coupled: if True, the weight is converted to (end-beg)*weight
    :type coupled: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16
    :param scaling: The scaling of the DFE
    :type scaling: float

    When coupled is True, the "weight" may be
    interpreted as a "per base pair" (or per unit, generally speaking) term.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    weight: float
    lo: float
    hi: float
    h: float = 1.0
    coupled: bool = True
    label: int = 0
    scaling: float = 1.0

    def __attrs_post_init__(self):
        super(UniformS, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.lo,
            self.hi,
            self.h,
            self.coupled,
            self.label,
            self.scaling,
        )


@_add_deprecated_properties
@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(eq=False, **_common_attr_attribs)
class MultivariateGaussianEffects(fwdpy11._fwdpy11._ll_MultivariateGaussianEffects):
    """
    Pleiotropic effects via a multivariate Gaussian distribution.

    This class can be used to generate mutations with both vectors
    of effect sizes as well as a separate fixed effect.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: Beginning of the region
    :type beg: float
    :param end: End of the region
    :type end: float
    :param weight: Weight on the region
    :type weight: float
    :param matrix: Variance-covariance matrix
    :type matrix: numpy.ndarray
    :param fixed_effect: Fixed effect size. Defaults to 0.0.
    :type fixed_effect: float
    :param h: Dominance. Defaults to 1.0
    :type h: float
    :param coupled: Specify if weight is function of end-beg or not. Defaults to True
    :type coupled: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16

    The input matrix must be square and semi-positive definite.   If either
    of these conditions are not met, ValueError will be raised. ValueError
    will also be raised if the input matrix contains any non-finite values.

    .. note::

        The dominance parameter (`h`) applies to both the fixed effect and those
        drawn from a multivariate normal.

    .. versionadded:: 0.3.0

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    beg: float
    end: float
    weight: float
    cov_matrix: np.ndarray
    fixed_effect: float = 0.0
    h: float = 1.0
    coupled: bool = True
    label: int = 0

    def __attrs_post_init__(self):
        super(MultivariateGaussianEffects, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.cov_matrix,
            self.fixed_effect,
            self.h,
            self.coupled,
            self.label,
        )

    def __eq__(self, other):
        if self.beg != other.beg:
            return False

        if self.beg != other.beg:
            return False
        if self.end != other.end:
            return False
        if self.weight != other.weight:
            return False
        if self.fixed_effect != other.fixed_effect:
            return False
        if self.h != other.h:
            return False
        if self.coupled != other.coupled:
            return False
        if self.label != other.label:
            return False
        return np.array_equal(self.cov_matrix, other.cov_matrix)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(eq=False, **_common_attr_attribs)
class mvDES(fwdpy11._fwdpy11._ll_mvDES):
    """
    General multivariate distribution of effect sizes.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param des: Distributions of effect sizes
    :type des: list
    :param means: means marginal gaussian Distributions
    :type means: numpy.ndarray
    :param matrix: Variance/covariance matrix
    :type matrix: numpy.ndarray

    .. versionadded:: 0.7.0

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    des: object
    means: object
    matrix: typing.Optional = None

    def __attrs_post_init__(self):
        super(mvDES, self).__init__(self.des, self.means, self.matrix)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(eq=False, **_common_attr_attribs)
class DiscreteDESD(fwdpy11._fwdpy11._ll_DiscreteDESD):
    """
    Discretized distribution of effect sizes and dominance.

    This class allows you to specify a discrete joint
    distrubtion of effect size and dominance.

    The distribution is specified by a list of tuples.
    Each tuple contains (effect size, dominance, weight).
    The weights must all be >= 0 and all values must fe finite.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param beg: Beginning of the region
    :type beg: float
    :param end: End of the region
    :type end: float
    :param weight: Weight on the region
    :type weight: float
    :param joint_dist: The joint distribution + weights
    :type joint_dist: list
    :param coupled: Specify if weight is function of end-beg or not. Defaults to True
    :type coupled: bool
    :param label: Fill :attr:`fwdpy11.Mutation.label` with this value.
    :type label: numpy.uint16
    :param scaling: The scaling of the DFE
    :type scaling: float

    .. versionadded:: 0.10.0
    """

    beg: float
    end: float
    weight: float
    joint_dist: typing.List[typing.Tuple[float, float, float]]
    coupled: bool = True
    label: int = 0
    scaling: float = 1.0

    def __attrs_post_init__(self):
        super(DiscreteDESD, self).__init__(
            self.beg,
            self.end,
            self.weight,
            self.joint_dist,
            self.coupled,
            self.label,
            self.scaling,
        )
