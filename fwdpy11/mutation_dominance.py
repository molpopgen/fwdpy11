#
# Copyright(C) 2021 Kevin Thornton < krthornt @uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software : you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.If not, see < http: //www.gnu.org/licenses/>.
#

import attr

import fwdpy11._fwdpy11

from .class_decorators import (attr_add_asblack, attr_class_pickle_with_super,
                               attr_class_to_from_dict)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class FixedDominance(fwdpy11._fwdpy11._ll_FixedDominance):
    """
    Fixed heterozygous effects.

    :param h: The heterozygous effect of a mutation, or "dominance".
    :type h: float

    .. versionadded:: 0.14.0
    """

    h: float = attr.ib(validator=attr.validators.instance_of(float))

    def __attrs_post_init__(self):
        super(FixedDominance, self).__init__(self.h)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class ExponentialDominance(fwdpy11._fwdpy11._ll_ExponentialDominance):
    """
    Exponential distribution of heterozygous effects.

    :param m: The mean of the distribution.
    :type m: float

    .. versionadded:: 0.14.0
    """

    m: float = attr.ib(validator=attr.validators.instance_of(float))

    def __attrs_post_init__(self):
        super(ExponentialDominance, self).__init__(self.m)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class UniformDominance(fwdpy11._fwdpy11._ll_UniformDominance):
    """
    Uniform distribution of heterozygous effects.

    :param lo: The lower bound of the range
    :type lo: float
    :param hi: The upper bound of the range
    :type hi: float

    .. versionadded:: 0.14.0
    """

    lo: float = attr.ib(validator=attr.validators.instance_of(float))
    hi: float = attr.ib(validator=attr.validators.instance_of(float))

    def __attrs_post_init__(self):
        super(UniformDominance, self).__init__(self.lo, self.hi)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11")
class LargeEffectExponentiallyRecessive(
    fwdpy11._fwdpy11._ll_LargeEffectExponentiallyRecessive
):
    """
    Large effect mutations are more recessive according to
    the function :math:`y \\times e^{-k|s|}`, where :math:`s` is the
    effect size and :math:`y` is a "scaling" parameter.

    :param k: The "tuning" parameter.
    :type k: float
    :type scaling: Scaling parameter
    :type scaling: float

    .. versionadded:: 0.14.0
    """

    k: float = attr.ib(validator=attr.validators.instance_of(float))
    scaling: float = attr.ib(validator=attr.validators.instance_of(float), default=1.0)

    @k.validator
    def _k_positive(self, attribute, value):
        if value <= 0.0:
            raise ValueError(f"{attribute} must be > 0.0")

    def __attrs_post_init__(self):
        super(LargeEffectExponentiallyRecessive, self).__init__(self.k, self.scaling)
