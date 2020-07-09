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
                               attr_class_to_from_dict,
                               attr_class_to_from_dict_no_recurse)

_common_attr_attribs = {"frozen": True, "auto_attribs": True, "repr_ns": "fwdpy11"}


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class Optimum(fwdpy11._fwdpy11._ll_Optimum):
    """
    Parameters for a trait optimum.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optimum: The trait value
    :type optimum: float
    :param VS: Strength of stabilizing selection
    :type VS: float
    :param when: The time when the optimum shifts
    :type when: int or None

    .. note::

        When used to model a stable optimum (e.g.,
        :class:`fwdpy11.GSS`), the ``when`` parameter is omitted.
        The ``when`` parameter is used for moving optima
        (:class:`fwdpy11.GSSmo`).

    .. versionadded:: 0.7.1

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    optimum: float = attr.ib(validator=attr.validators.instance_of(float))
    VS: float = attr.ib(validator=attr.validators.instance_of(float))
    when: typing.Optional[int] = attr.ib(default=None)

    @when.validator
    def validate_when(self, attribute, value):
        if value is not None:
            attr.validators.instance_of(int)(self, attribute, value)

    def __attrs_post_init__(self):
        super(Optimum, self).__init__(self.optimum, self.VS, self.when)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs, eq=False)
class PleiotropicOptima(fwdpy11._fwdpy11._ll_PleiotropicOptima):
    """
    Parameters for multiple trait optima

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optima: The trait values
    :type optima: List[float]
    :param VS: Strength of stabilizing selection
    :type VS: float
    :param when: The time when the optimum shifts
    :type when: int or None

    .. note::

        When used to model stable optima (e.g.,
        :class:`fwdpy11.MultivariateGSS`), the ``when`` parameter is omitted.
        The ``when`` parameter is used for moving optima
        (:class:`fwdpy11.MultivariateGSSmo`).

    .. versionadded:: 0.7.1

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    optima: typing.List[float]
    VS: float = attr.ib(validator=attr.validators.instance_of(float))
    when: typing.Optional[int] = attr.ib(default=None)

    @when.validator
    def validate_when(self, attribute, value):
        if value is not None:
            attr.validators.instance_of(int)(self, attribute, value)

    def __attrs_post_init__(self):
        super(PleiotropicOptima, self).__init__(self.optima, self.VS, self.when)

    def __eq__(self, other):
        optima_equal = np.array_equal(self.optima, other.optima)
        VS_equal = self.VS == other.VS
        when_equal = False
        if self.when is not None and other.when is not None:
            when_equal = self.when == other.when

        return optima_equal and VS_equal and when_equal


@attr_add_asblack
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class GSS(fwdpy11._fwdpy11._ll_GSSmo):
    """
    Gaussian stabilizing selection on a single trait.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optimum: The optimal trait value
    :type optimum: float or fwdpy11.Optimum
    :param VS: Inverse strength of stabilizing selection
    :type VS: float or None

    .. note::

        VS should be None if optimum is an instance
        of :class:`fwdpy11.Optimum`

    .. versionchanged:: 0.7.1

        Allow instances of fwdpy11.Optimum for intitialization

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    optimum: typing.Union[Optimum, float]
    VS: typing.Optional[float] = None

    def __attrs_post_init__(self):
        if self.VS is None:
            super(GSS, self).__init__(
                [Optimum(optimum=self.optimum.optimum, VS=self.optimum.VS, when=0)]
            )
        else:
            super(GSS, self).__init__(
                [Optimum(optimum=self.optimum, VS=self.VS, when=0)]
            )

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        if self.VS is None:
            super(GSS, self).__init__(
                [Optimum(optimum=self.optimum.optimum, VS=self.optimum.VS, when=0)]
            )
        else:
            super(GSS, self).__init__(
                [Optimum(optimum=self.optimum, VS=self.VS, when=0)]
            )


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class GSSmo(fwdpy11._fwdpy11._ll_GSSmo):
    """
    Gaussian stabilizing selection on a single trait with moving
    optimum.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optima: The optimal trait values
    :type optima: list[fwdpy11.Optimum]

    .. note::

        Instances of fwdpy11.Optimum must have valid
        values for ``when``.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class

    """

    optima: typing.List[Optimum] = attr.ib()

    @optima.validator
    def validate_optima(self, attribute, value):
        if len(value) == 0:
            raise ValueError("list of optima cannot be empty")
        for o in value:
            if o.when is None:
                raise ValueError("Optimum.when is None")

    def __attrs_post_init__(self):
        super(GSSmo, self).__init__(self.optima)


@attr_add_asblack
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class MultivariateGSS(fwdpy11._fwdpy11._ll_MultivariateGSSmo):
    """
    Multivariate gaussian stablizing selection.

    Maps a multidimensional trait to fitness using the Euclidian
    distance of a vector of trait values to a vector of optima.

    Essentially, this is Equation 1 of

    Simons, Yuval B., Kevin Bullaughey, Richard R. Hudson, and Guy Sella. 2018.
    "A Population Genetic Interpretation of GWAS Findings for Human Quantitative Traits."
    PLoS Biology 16 (3): e2002985.

    For the case of moving optima, see :class:`fwdpy11.MultivariateGSSmo`.
    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optima: The optimum value for each trait over time
    :type optima: numpy.array or list[fwdpy11.PleiotropicOptima]
    :param VS: Inverse strength of stablizing selection
    :type VS: float or None

    .. note::

        ``VS`` should be ``None`` if ``optima`` is list[fwdpy11.PleiotropicOptima]

        ``VS`` is :math:`\omega^2` in the Simons et al. notation

    .. versionchanged:: 0.7.1
        Allow initialization with list of fwdpy11.PleiotropicOptima

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    optima: typing.Union[PleiotropicOptima, typing.List[float]]
    VS: typing.Optional[float] = None

    def __attrs_post_init__(self):
        if self.VS is None:
            super(MultivariateGSS, self).__init__([self.optima])
        else:
            super(MultivariateGSS, self).__init__(self._convert_to_list())

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
        if self.VS is None:
            super(MultivariateGSS, self).__init__([self.optima])
        else:
            super(MultivariateGSS, self).__init__(self._convert_to_list())

    def _convert_to_list(self):
        if self.VS is None:
            raise ValueError("VS must not be None")
        return [PleiotropicOptima(optima=self.optima, VS=self.VS, when=0)]


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class MultivariateGSSmo(fwdpy11._fwdpy11._ll_MultivariateGSSmo):
    """
    Multivariate gaussian stablizing selection with moving optima
    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optima: list of optima over time
    :type optima: list[fwdpy11.PleiotropicOptima]

    .. versionchanged:: 0.7.1

        Allow initialization with list of fwdpy11.PleiotropicOptima

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    optima: typing.List[PleiotropicOptima] = attr.ib()

    @optima.validator
    def validate_optima(self, attribute, value):
        if len(value) == 0:
            raise ValueError("list of optima cannot be empty")
        for o in value:
            if o.when is None:
                raise ValueError("PleiotropicOptima.when is None")

    def __attrs_post_init__(self):
        super(MultivariateGSSmo, self).__init__(self.optima)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class NoNoise(fwdpy11._fwdpy11._ll_NoNoise):
    """
    No random effects on genetic values

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    def __attrs_post_init__(self):
        super(NoNoise, self).__init__()


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class GaussianNoise(fwdpy11._fwdpy11._ll_GaussianNoise):
    """
    Gaussian noise added to genetic values.
    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param sd: Standard deviation
    :type sd: float
    :param mean: Mean value
    :type mean: float

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    sd: float
    mean: float = 0.0

    def __attrs_post_init__(self):
        super(GaussianNoise, self).__init__(self.sd, self.mean)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class Additive(fwdpy11._fwdpy11._ll_Additive):
    """
    Additive effects on genetic values.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param scaling: How to treat mutant homozygotes.
    :type scaling: float
    :param gvalue_to_fitness: How to map trait value to fitness
    :type gvalue_to_fitness: fwdpy11.GeneticValueIsTrait
    :param noise: Random effects on trait values
    :type noise: fwdpy11.GeneticValueNoise

    When ``gvalue_to_fitness`` is ``None``, then we are
    modeling additive effects on fitness.

    For a model of fitness, the genetic value is 1, 1+e*h,
    1+``scaling``*e for genotypes AA, Aa, and aa, respectively,
    where ``e`` and ``h`` are the effect size and dominance, respectively.

    For a model of a trait (phenotype), meaning ``gvalue_to_fitness``
    is not ``None``, the values for the three genotypes are 0, e*h,
    and e, respectively.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    scaling: float
    gvalue_to_fitness: fwdpy11.GeneticValueIsTrait = None
    noise: fwdpy11.GeneticValueNoise = None
    ndemes: int = 1

    def __attrs_post_init__(self):
        super(Additive, self).__init__(
            self.scaling, self.gvalue_to_fitness, self.noise, self.ndemes
        )


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class Multiplicative(fwdpy11._fwdpy11._ll_Multiplicative):
    """
    Multiplicative effects on genetic values.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param scaling: How to treat mutant homozygotes.
    :type scaling: float
    :param gvalue_to_fitness: How to map trait value to fitness
    :type gvalue_to_fitness: fwdpy11.GeneticValueIsTrait
    :param noise: Random effects on trait values
    :type noise: fwdpy11.GeneticValueNoise

    When ``gvalue_to_fitness`` is ``None``, then we are
    modeling multiplicative effects on fitness.

    For a model of fitness, the genetic value is 1, 1+e*h,
    1+``scaling``*e for genotypes AA, Aa, and aa, respectively,
    where ``e`` and ``h`` are the effect size and dominance, respectively.

    For a model of a trait (phenotype), meaning ``gvalue_to_fitness``
    is not ``None``, the values for the three genotypes are 0, e*h,
    and e, respectively.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    scaling: float
    gvalue_to_fitness: fwdpy11.GeneticValueIsTrait = None
    noise: fwdpy11.GeneticValueNoise = None
    ndemes: int = 1

    def __attrs_post_init__(self):
        super(Multiplicative, self).__init__(
            self.scaling, self.gvalue_to_fitness, self.noise, self.ndemes
        )


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class GBR(fwdpy11._fwdpy11._ll_GBR):
    """
    The "gene-based recessive" trait model described in Thornton et al.
    2013 http://dx.doi.org/10.1371/journal.pgen.1003258 and Sanjak et al. 2017
    http://dx.doi.org/10.1371/journal.pgen.1006573.

    The trait value is the geometric mean of the sum of effect sizes on each haplotype.
    It is undefined for the case where these sums are negative.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param gvalue_to_fitness: How to map trait value to fitness
    :type gvalue_to_fitness: fwdpy11.GeneticValueIsTrait
    :param noise: Random effects on trait values
    :type noise: fwdpy11.GeneticValueNoise

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    gvalue_to_fitness: object
    noise: object = None

    def __attrs_post_init__(self):
        super(GBR, self).__init__(self.gvalue_to_fitness, self.noise)


@attr_add_asblack
@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class StrictAdditiveMultivariateEffects(
    fwdpy11._fwdpy11._ll_StrictAdditiveMultivariateEffects
):
    """
    Multivariate trait values under strictly additive effects.

    Calculate the trait value for a diploid in a :class:`fwdpy11.DiploidPopulation`
    for a multidimensional trait.

    This class is restricted to the case of simple additive effects, meaning
    that any dominance terms associated with mutations are ignored.

    During a simulation, :attr:`fwdpy11.DiploidMetadata.g` is filled with the
    genetic value corresponding to a "focal" trait specified upon object construction.
    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param ndimensions: Number of trait dimensions
    :type ndimensions: int
    :param focal_trait: Index of the focal trait
    :type focal_trait: int
    :param gvalue_to_fitness: Function mapping trait value to fitness
    :type gvalue_to_fitness: :class:`fwdpy11.GeneticValueToFitnessMap`
    :param noise: Function adding random additive noise to trait value
    :type noise: :class:`fwdpy11.GeneticValueNoise`

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    ndimensions: int
    focal_trait: int
    gvalue_to_fitness: object
    noise: object = None

    def __attrs_post_init__(self):
        super(StrictAdditiveMultivariateEffects, self).__init__(
            self.ndimensions, self.focal_trait, self.gvalue_to_fitness, self.noise
        )
