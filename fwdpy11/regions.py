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
import math
import numpy as np


class Region(object):
    """
    Representation of a "region" in a simulation.

    Attributes:
        b: the beginning of the region
        e: the end of the region
        w: the "weight" assigned to the region
        l: A label assigned to the region.
            Labels must be integers, and can be used to
            'tag' mutations arising in different regions.

    See :func:`evolve_regions` for how this class may be used to
        parameterize a simulation.

    This class is extended by:
        * :class:`fwdpy11.Sregion`
    """

    def __init__(self, beg, end, weight, coupled=True, label=0):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used to take mutations
            from this region.

        When coupled is True, the "weight" may be interpreted
        as a "per base pair"
        (or per unit, generally speaking) term.

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            r = fwdpy11.Region(0,1,1)
            #A more "biological" case:
            #  The region covers positions 1 through 1,000,
            #  and the per-base pair "weight" is 1e-5:
            r = fwdpy11.Region(1,1000,1e-5,True)
        """
        if math.isinf(beg):
            raise ValueError("fwdpy11.Region: beg not finite")
        if math.isinf(end):
            raise ValueError("fwdpy11.Region: end not finite")
        if math.isinf(weight):
            raise ValueError("fwdpy11.Region: weight not finite")
        if math.isnan(beg):
            raise ValueError("fwdpy11.Region: beg not a number")
        if math.isnan(end):
            raise ValueError("fwdpy11.Region: end not a number")
        if math.isnan(weight):
            raise ValueError("fwdpy11.Region: weight not a number")
        if weight < 0.0:
            raise ValueError("fwdpy11.Region: weight < 0.0")
        self.b = float(beg)
        self.e = float(end)
        self.w = float(weight)
        self.c = coupled
        self.l = label
        if self.c is True:
            self.w = (self.e - self.b) * self.w

    def __repr__(self):
        x = 'regions.Region(beg=%s,end=%s,'
        x += 'weight=%s,coupled=%s,label=%s)'
        return x % (self.b, self.e, self.w, self.c, self.l)


class Sregion(Region):
    """
    Representation of a "region" in a simulation with a dominance term.

    This class is the base class for a general set
    of objects representing distributions of fitness effects.

    Attributes:
        * b: the beginning of the region
        * e: the end of the region
        * w: the "weight" assigned to the region
        * h: the dominance term
        * l: A label assigned to the region.
            Labels must be integers, and can be used to
            'tag' mutations arising in different regions.
        * scaling: The scaling of the distrubution.  See note below.

    See :func:`evolve_regions` for how this class
    may be used to parameterize a simulation.

    .. note:: This class cannot be used directly to parameterize a simulation.
        Rather, you must use a derived type that specifies a
        distribution of fitness effects.  These types include:
       :class:`fwdpy11.fwdpy11.ConstantS`,
       :class:`fwdpy11.fwdpy11.UniformS`,
       :class:`fwdpy11.fwdpy11.ExpS`,
       :class:`fwdpy11.fwdpy11.GammaS`, and
       :class:`fwdpy11.fwdpy11.GaussianS`

    .. note:: The scaling of a distribution refers to the distribution of effect sizes.
        For example, if scaling = 1.0, then the DFE is on the effect size itself.  If
        scaling = 2N (where N is the population size), then the DFE is on 2Ns.  If N
        is not constant during a simulation, then the scaling is with respect to some
        "reference" population size.

    .. versionchanged:: 0.13.a2
        Added "scaling" attribute.

    """

    def __init__(self, beg, end, weight, h=1.0, coupled=True, label=0, scaling=1):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param h: the dominance
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used to take mutations
            from this region.
        :param scaling: The scaling of the DFE

        When coupled is True, the "weight" may be
        interpreted as a "per base pair"
        (or per unit, generally speaking) term.

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            #Examples for models where the 3 genotype fitnesses are
            #1, 1+sh, and 1+2s, respectively
            recessive = fwdpy11.Sregion(0,1,1,0)
            additive = fwdpy11.Sregion(0,1,1,1.0)
            dominant = fwdpy11.Sregion(0,1,1,2.0)
        """
        if math.isinf(h):
            raise ValueError("fwdpy11.Sregion: h not finite")
        if math.isnan(h):
            raise ValueError("fwdpy11.Segion: h not a number")
        self.h = float(h)
        self.scaling = np.uint32(scaling)
        super(Sregion, self).__init__(beg, end, weight, coupled, label)

    # def __str__(self):
    #    retval = "h = " + "{:.9f}".format(self.h)
    #    retval += ", " + super(Sregion, self).__str__()


class GammaS(Sregion):
    """
    Gamma distribution of fitness effects

    Attributes:
        * b: the beginning of the region
        * e: the end of the region
        * w: the "weight" assigned to the region
        * mean: mean of the Gamma
        * shape: shape of the Gamma
        * h: the dominance term
        * l: A label assigned to the region.
            Labels must be integers, and can be used to
            'tag' mutations arising in different regions.
        * scaling: The scaling of the distrubution.

    See :func:`evolve_regions` for how this
    class may be used to parameterize a simulation
    """

    def __init__(self, beg, end, weight, mean, shape, h=1.0,
                 coupled=True, label=0, scaling=1):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param mean: the mean selection coefficient
        :param shape: the shape parameter of the distribution
        :param h: the dominance
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used
            to take mutations from this region.
        :param scaling: The scaling of the DFE

        When coupled is True, the "weight" may be
        interpreted as a "per base pair"
        (or per unit, generally speaking) term.

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            gdist = fwdpy11.GammaS(0,1,1,-0.1,0.35)
        """
        if math.isinf(mean):
            raise ValueError("fwdpy11.GammaS: mean not finite")
        if math.isinf(shape):
            raise ValueError("fwdpy11.GammaS: shape not finite")
        if math.isnan(mean):
            raise ValueError("fwdpy11.GammaS: mean not a number")
        if math.isnan(shape):
            raise ValueError("fwdpy11.GammaS: shape not a number")
        self.mean = float(mean)
        self.shape = float(shape)
        super(GammaS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)

    def callback(self):
        """
        Returns the C++ callback for this DFE
        """
        from .fwdpp_extensions import makeGammaSH
        return makeGammaSH(self.mean, self.shape, self.h, self.scaling)

    def __repr__(self):
        x = 'regions.GammaS(beg=%s,end=%s,weight=%s,'
        x += 'mean=%s,shape=%s,coupled=%s,label=%s,scaling=%s)'
        return x % (self.b, self.e, self.w, self.mean,
                    self.shape, self.c, self.l, self.scaling)


class ConstantS(Sregion):
    """
    Constant/fixed selection coefficient

    Attributes:
        * b: the beginning of the region
        * e: the end of the region
        * w: the "weight" assigned to the region
        * s: the selection coefficient
        * h: the dominance term
        * l: A label assigned to the region.
            Labels must be integers, and can be used to
            'tag' mutations arising in different regions.
        * scaling: The scaling of the distrubution.

    See :func:`evolve_regions` for how this class may be used to
    parameterize a simulation
    """

    def __init__(self, beg, end, weight, s, h=1.0, coupled=True, label=0, scaling=1):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param s: the selection coefficient
        :param h: the dominance
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used to take mutations
            from this region.
        :param scaling: The scaling of the DFE

        When coupled is True, the "weight" may be interpreted
        as a "per base pair"
        (or per unit, generally speaking) term.

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            #s = -0.1 and h = 0
            constantS = fwdpy11.ConstantS(0,1,1,-0.1,0)
        """
        if math.isinf(s):
            raise ValueError("fwdpy11.ConstantS: s not finite")
        if math.isnan(s):
            raise ValueError("fwdpy11.ConstantS: s not a number")
        self.s = float(s)
        super(ConstantS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)

    def callback(self):
        """
        Returns the C++ callback for this DFE
        """
        from .fwdpp_extensions import makeConstantSH
        return makeConstantSH(self.s, self.h, self.scaling)

    def __repr__(self):
        x = 'regions.ConstantS(beg=%s,end=%s,weight=%s,'
        x += 's=%s,h=%s,coupled=%s,label=%s,scaling=%s)'
        return x % (self.b, self.e, self.w, self.s, self.h, self.c, self.l, self.scaling)


class UniformS(Sregion):
    """
    Uniform distribution on selection coefficients

    Attributes:
        * b: the beginning of the region
        * e: the end of the region
        * w: the "weight" assigned to the region
        * lo: the lower bound on s
        * hi: the upper bound on s
        * h: the dominance term
        * l: A label assigned to the region.
            Labels must be integers, and can be used to
            'tag' mutations arising in different regions.
        * scaling: The scaling of the distrubution.

    See :func:`evolve_regions` for how this
    class may be used to parameterize a simulation
    """

    def __init__(self, beg, end, weight, lo, hi, h=1.0, coupled=True, label=0, scaling=1):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param lo: lower bound on s
        :param hi: upper bound on s
        :param h: the dominance
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used to take
            mutations from this region.
        :param scaling: The scaling of the DFE

        When coupled is True, the "weight" may be
        interpreted as a "per base pair"
        (or per unit, generally speaking) term.

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            #s is uniform on [0,-1]
            constantS = fwdpy11.UniformS(0,1,1,0,-1,0)
        """
        if math.isinf(lo):
            raise ValueError("fwdpy11.UniformS: lo not finite")
        if math.isinf(hi):
            raise ValueError("fwdpy11.UniformS: hi not finite")
        if math.isnan(lo):
            raise ValueError("fwdpy11.UniformS: lo not a number")
        if math.isnan(hi):
            raise ValueError("fwdpy11.UniformS: hi not a number")
        self.lo = float(lo)
        self.hi = float(hi)
        super(UniformS, self).__init__(beg, end, weight, h, coupled, scaling)

    def callback(self):
        """
        Returns the C++ callback for this DFE
        """
        from .fwdpp_extensions import makeUniformSH
        return makeUniformSH(self.lo, self.hi, self.h, self.scaling)

    def __repr__(self):
        x = 'regions.UniformS(beg=%s,end=%s,weight=%s,'
        x += 'lo=%s,hi=%s,h=%s,coupled=%s,label=%s,scaling=%s)'
        return x % (self.b, self.e, self.w, self.lo,
                    self.hi, self.h, self.c, self.l, self.scaling)


class ExpS(Sregion):
    """
    Exponential distribution on selection coefficients

    Attributes:
        * b: the beginning of the region
        * e: the end of the region
        * w: the "weight" assigned to the region
        * mean: the mean selection coefficient
        * h: the dominance term
        * l: A label assigned to the region.
            Labels must be integers, and can be used to
            'tag' mutations arising in different regions.
        * scaling: The scaling of the distrubution.

    See :func:`evolve_regions` for how this class may be used to
    parameterize a simulation
    """

    def __init__(self, beg, end, weight, mean, h=1.0, coupled=True, label=0, scaling=1):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param mean: the mean selection coefficient
        :param h: the dominance
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used to take
            mutations from this region.
        :param scaling: The scaling of the DFE

        When coupled is True, the "weight" may be
        interpreted as a "per base pair"
        (or per unit, generally speaking) term.

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            #s is exp(-0.1) and recessive
            constantS = fwdpy11.ExpS(0,1,1,0,-0.1,0)
        """
        if math.isinf(mean):
            raise ValueError("fwdpy11.ExpS: mean not finite")
        if math.isnan(mean):
            raise ValueError("fwdpy11.ExpS: mean not a number")
        self.mean = float(mean)
        super(ExpS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)

    def callback(self):
        """
        Returns the C++ callback for this DFE
        """
        from .fwdpp_extensions import makeExpSH
        return makeExpSH(self.mean, self.h, self.scaling)

    def __repr__(self):
        x = 'regions.ExpS(beg=%s,end=%s,weight=%s,'
        x += 'mean=%s,h=%s,coupled=%s,label=%s,scaling=%s)'
        return x % (self.b, self.e, self.w, self.mean, self.h, self.c, self.l, self.scaling)


class GaussianS(Sregion):
    """
    Gaussian distribution on selection coefficients
    (effect sizes for sims of quantitative traits)

    Attributes:
        * b: the beginning of the region
        * e: the end of the region
        * w: the "weight" assigned to the region
        * sd: the standard deviation
        * h: the dominance ter
        * l: A label assigned to the region.
            Labels must be integers, and can be used to
            'tag' mutations arising in different regions.
        * scaling: The scaling of the distrubution.

    The mean is zero.

    See :func:`evolve_regions` for how this class may be used to
    parameterize a simulation
    """

    def __init__(self, beg, end, weight, sd, h=1.0, coupled=True, label=0, scaling=1):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param mean: the mean selection coefficient
        :param h: the dominance
        :param coupled: if True, the weight is converted to (end-beg)*weight
        :param label: Not relevant to recombining regions.
            Otherwise, this value will be used to take mutations
            from this region.
        :param scaling: The scaling of the DFE

        When coupled is True, the "weight" may be
        interpreted as a "per base pair"
        (or per unit, generally speaking) term.

        Example:

        .. testcode::

            #A simple case
            import fwdpy11
            #s N(0,0.1) and co-dominant
            constantS = fwdpy11.GaussianS(0,1,1,0,0.1,1)
        """
        if math.isinf(sd):
            raise ValueError("fwdpy11.GaussianS: sd not finite")
        if math.isnan(sd):
            raise ValueError("fwdpy11.GaussianS: sd not a number")
        self.sd = float(sd)
        super(GaussianS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)

    def callback(self):
        """
        Returns the C++ callback for this DFE
        """
        from .fwdpp_extensions import makeGaussianSH
        return makeGaussianSH(self.sd, self.h, self.scaling)

    def __repr__(self):
        x = 'regions.GaussianS(beg=%s,end=%s,weight=%s,'
        x += 'sd=%s,h=%s,coupled=%s,label=%s,scaling=%s)'
        return x % (self.b, self.e, self.w, self.sd, self.h, self.c, self.l, self.scaling)
