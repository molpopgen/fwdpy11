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

import enum
import typing
import warnings

import attr
import demes
import numpy as np


from ._fwdpy11 import (
    GeneticValueIsTrait,
    GeneticValueNoise,
    _ll_Additive,
    _ll_GaussianNoise,
    _ll_GBR,
    _ll_GaussianStabilizingSelection,
    _ll_GSSmo,
    _ll_Multiplicative,
    _ll_MultivariateGSSmo,
    _ll_NoNoise,
    _ll_Optimum,
    _ll_PleiotropicOptima,
    _ll_StrictAdditiveMultivariateEffects,
    _PyDiploidGeneticValue,
)
from .class_decorators import (
    attr_class_pickle_with_super,
    attr_class_to_from_dict,
    attr_class_to_from_dict_no_recurse,
)

from ._types.forward_demes_graph import ForwardDemesGraph


class TimingError(Exception):
    pass


@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, auto_detect=True)
class Optimum(_ll_Optimum):
    """
    Parameters for a trait optimum.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optimum: The trait value
    :type optimum: float
    :param VS: Strength of stabilizing selection
    :type VS: float
    :param when: The time when the optimum shifts
    :type when: int or None

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

    def __repr__(self) -> str:
        return (
            f"fwdpy11.Optimum(optimum={self.optimum}, VS={self.VS}, when={self.when})"
        )


@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True, eq=False, auto_detect=True)
class PleiotropicOptima(_ll_PleiotropicOptima):
    """
    Parameters for multiple trait optima

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param optima: The trait values
    :type optima: List[float]
    :param VS: Strength of stabilizing selection
    :type VS: float
    :param when: The time when the optimum shifts
    :type when: int or None

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

    def __repr__(self) -> str:
        args = f"optima={self.optima}, VS={self.VS}, when={self.when}"
        return f"(fwdpy11.PleiotropicOptima({args}))"


@attr.s(auto_attribs=True, frozen=True)
class GaussianStabilizingSelection(_ll_GaussianStabilizingSelection):
    """
    Define a mapping of phenotype-to-fitness according to a
    Gaussian stabilizing selection model.

    Instances of this trait must be constructed by one of the
    various class methods available.
    """

    is_single_trait: bool
    optima: list

    def __attrs_post_init__(self):
        if self.optima != sorted(self.optima, key=lambda x: x.when):
            raise ValueError("optima must be sorted by time from past to present")
        if self.is_single_trait is True:
            inner = _ll_GSSmo(self.optima)
        else:
            inner = _ll_MultivariateGSSmo(self.optima)
        super(GaussianStabilizingSelection, self).__init__(inner)

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(**d)
        self.__attrs_post_init__()

    @classmethod
    def single_trait(cls, optima: typing.List[Optimum]):
        """
        Stabilizing selection on a single trait

        :param optima: The optimum values.
                       Multiple values specify a moving optimum.
        :type optima: List[fwdpy11.Optimum]

        """
        new_optima = []
        for o in optima:
            if o.when is None:
                new_optima.append(Optimum(o.optimum, o.VS, 0))
            else:
                new_optima.append(o)
        return cls(is_single_trait=True, optima=new_optima)

    @classmethod
    def pleiotropy(cls, optima: typing.List[PleiotropicOptima]):
        """
        Stabilizing selection with pleiotropy.

        :param optima: The optimum values.
                       Multiple values specify a moving optimum.
        :type optima: List[fwdpy11.PleiotropicOptima]
        """
        return cls(is_single_trait=False, optima=optima)

    def asdict(self):
        return {"is_single_trait": self.is_single_trait, "optima": self.optima}

    @classmethod
    def fromdict(cls, data):
        return cls(**data)

    def validate_timings(self, deme: int, demography: ForwardDemesGraph) -> None:
        graph = demography.demes_graph
        deme_obj = graph.demes[deme]
        if deme_obj.start_time == float("inf"):
            start_time = demography.to_backwards_time(0)
        else:
            start_time = deme_obj.start_time
        for i in self.optima[:1]:
            btime = demography.to_backwards_time(i.when)
            if btime < start_time:
                msg = f"deme {deme_obj.name}, event {i} is at time {btime} in the past "
                msg += f", but the deme starts at {int(start_time)}"
                raise TimingError(msg)


@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True)
class NoNoise(_ll_NoNoise):
    """
    No random effects on genetic values

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    def __attrs_post_init__(self):
        super(NoNoise, self).__init__()


@attr_class_pickle_with_super
@attr_class_to_from_dict
@attr.s(auto_attribs=True, frozen=True)
class GaussianNoise(_ll_GaussianNoise):
    """
    Gaussian noise added to genetic values.
    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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


@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(auto_attribs=True, frozen=True)
class Additive(_ll_Additive):
    """
    Additive effects on genetic values.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param scaling: How to treat mutant homozygotes.
    :type scaling: float
    :param gvalue_to_fitness: How to map trait value to fitness
    :type gvalue_to_fitness: fwdpy11.GeneticValueIsTrait
    :param noise: Random effects on trait values
    :type noise: fwdpy11.GeneticValueNoise

    When `gvalue_to_fitness` is `None`, then we are
    modeling additive effects on fitness.

    For a model of fitness, the genetic value is 1, 1+e*h,
    1+`scaling`*e for genotypes AA, Aa, and aa, respectively,
    where `e` and `h` are the effect size and dominance, respectively.

    For a model of a trait (phenotype), meaning `gvalue_to_fitness`
    is not `None`, the values for the three genotypes are 0, e*h,
    and e, respectively.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    scaling: float
    gvalue_to_fitness: GeneticValueIsTrait = None
    noise: GeneticValueNoise = None
    ndemes: int = 1

    def __attrs_post_init__(self):
        super(Additive, self).__init__(
            self.scaling, self.gvalue_to_fitness, self.noise, self.ndemes
        )

    def validate_timings(self, deme: int, demography: ForwardDemesGraph) -> None:
        if self.gvalue_to_fitness is None:
            return
        self.gvalue_to_fitness.validate_timings(deme, demography)


@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(auto_attribs=True, frozen=True)
class Multiplicative(_ll_Multiplicative):
    """
    Multiplicative effects on genetic values.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param scaling: How to treat mutant homozygotes.
    :type scaling: float
    :param gvalue_to_fitness: How to map trait value to fitness
    :type gvalue_to_fitness: fwdpy11.GeneticValueIsTrait
    :param noise: Random effects on trait values
    :type noise: fwdpy11.GeneticValueNoise

    When `gvalue_to_fitness` is `None`, then we are
    modeling multiplicative effects on fitness.

    For a model of fitness, the genetic value is 1, 1+e*h,
    1+`scaling`*e for genotypes AA, Aa, and aa, respectively,
    where `e` and `h` are the effect size and dominance, respectively.

    For a model of a trait (phenotype), meaning `gvalue_to_fitness`
    is not `None`, the values for the three genotypes are 0, e*h,
    and e, respectively.

    .. versionchanged:: 0.8.0

        Refactored to use attrs and inherit from
        low-level C++ class
    """

    scaling: float
    gvalue_to_fitness: GeneticValueIsTrait = None
    noise: GeneticValueNoise = None
    ndemes: int = 1

    def __attrs_post_init__(self):
        super(Multiplicative, self).__init__(
            self.scaling, self.gvalue_to_fitness, self.noise, self.ndemes
        )

    def validate_timings(self, deme: int, demography: ForwardDemesGraph) -> None:
        if self.gvalue_to_fitness is None:
            return
        self.gvalue_to_fitness.validate_timings(deme, demography)


@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(auto_attribs=True, frozen=True)
class GBR(_ll_GBR):
    """
    The "gene-based recessive" trait model described in Thornton et al.
    2013 http://dx.doi.org/10.1371/journal.pgen.1003258 and Sanjak et al. 2017
    http://dx.doi.org/10.1371/journal.pgen.1006573.

    The trait value is the geometric mean of the sum of effect
    sizes on each haplotype.
    It is undefined for the case where these sums are negative.

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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

    def validate_timings(self, deme: int, demography: ForwardDemesGraph) -> None:
        if self.gvalue_to_fitness is None:
            return
        self.gvalue_to_fitness.validate_timings(deme, demography)


@attr_class_pickle_with_super
@attr_class_to_from_dict_no_recurse
@attr.s(auto_attribs=True, frozen=True)
class AdditivePleiotropy(_ll_StrictAdditiveMultivariateEffects):
    """
    Multivariate trait values under strictly additive effects.

    Calculate the trait value for a diploid in a
    :class:`fwdpy11.DiploidPopulation` for a multidimensional trait.

    This class is restricted to the case of simple additive effects, meaning
    that any dominance terms associated with mutations are ignored.

    During a simulation, :attr:`fwdpy11.DiploidMetadata.g` is filled with the
    genetic value corresponding to a "focal" trait specified upon
    object construction. This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
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
        super(AdditivePleiotropy, self).__init__(
            self.ndimensions, self.focal_trait, self.gvalue_to_fitness, self.noise
        )

    def validate_timings(self, deme: int, demography: ForwardDemesGraph) -> None:
        if self.gvalue_to_fitness is None:
            return
        self.gvalue_to_fitness.validate_timings(deme, demography)


class PyDiploidGeneticValue(_PyDiploidGeneticValue):
    def __init__(self, ndim: int, gvalue_to_fitness=None, noise=None):
        self.gvalue_to_fitness = gvalue_to_fitness
        self.noise = noise
        super(PyDiploidGeneticValue, self).__init__(
            ndim,
            gvalue_to_fitness,
            noise,
        )

    def validate_timings(self, deme: int, demography: ForwardDemesGraph) -> None:
        if self.gvalue_to_fitness is None:
            return
        self.gvalue_to_fitness.validate_timings(deme, demography)
