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


class Region(object):
    """
    Representation of a "region" in a simulation.

    Attributes:
        * b: the beginning of the region
        * e: the end of the region
        * w: the "weight" assigned to the region
        * l: A label assigned to the region.  Labels must be integers,
            and can be used to 'tag' mutations arising in different regions.

    See :func:`fwdpy11.wright_fisher.evolve` for how this class may be used to
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
        self.b = float(beg)
        self.e = float(end)
        self.w = float(weight)
        self.c = coupled
        self.l = label
        if self.c is True:
            self.w = (self.e - self.b) * self.w

    @property
    def b(self):
        """
        Beginning of a region
        """
        return self.__b

    @b.setter
    def b(self, beg):
        if math.isinf(beg):
            raise ValueError("fwdpy11.Region: beg not finite")
        if math.isnan(beg):
            raise ValueError("fwdpy11.Region: beg not a number")
        self.__b = beg

    @property
    def e(self):
        """
        End of a region
        """
        return self.__e

    @e.setter
    def e(self, end):
        if math.isnan(end):
            raise ValueError("fwdpy11.Region: end not a number")
        if math.isinf(end):
            raise ValueError("fwdpy11.Region: end not finite")
        self.__e = end

    @property
    def w(self):
        """
        Weight on a region.
        """
        return self.__w

    @w.setter
    def w(self, weight):
        if math.isinf(weight):
            raise ValueError("fwdpy11.Region: weight not finite")
        if math.isnan(weight):
            raise ValueError("fwdpy11.Region: weight not a number")
        if weight < 0.0:
            raise ValueError("fwdpy11.Region: weight < 0.0")
        self.__w = weight

    @property
    def c(self):
        """
        If weight is coupled to end-beg. (boolean)
        """
        return self.__c

    @c.setter
    def c(self, coupled):
        self.__c = coupled

    @property
    def l(self):
        """
        Region label.
        """
        return self.__l

    @l.setter
    def l(self, label):
        import numpy as np
        self.__l = np.uint16(label)

    def __repr__(self):
        x = 'Region(beg=%s,end=%s,'
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
            Labels must be integers, and can be
            used to 'tag' mutations arising in different regions.
        * scaling: The scaling of the distrubution.  See note below.

    See :func:`fwdpy11.wright_fisher.evolve` for how this class
    may be used to parameterize a simulation.

    .. note:: This class cannot be used directly to parameterize a simulation.
        Rather, you must use a derived type that specifies a
        distribution of fitness effects.  These types include:
       :class:`fwdpy11.ConstantS`,
       :class:`fwdpy11.UniformS`,
       :class:`fwdpy11.ExpS`,
       :class:`fwdpy11.GammaS`, and
       :class:`fwdpy11.GaussianS`

    .. note:: The scaling of a distribution refers to the distribution of effect sizes.
        For example, if scaling = 1.0, then the DFE is on the effect size itself.  If
        scaling = 2N (where N is the population size), then the DFE is on 2Ns.  If N
        is not constant during a simulation, then the scaling is with respect to some
        "reference" population size.

    .. versionchanged:: 0.13.a2
        Added "scaling" attribute.

    """

    def __init__(self, beg, end, weight, h=1.0,
                 coupled=True, label=0, scaling=1):
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
        super(Sregion, self).__init__(beg, end, weight, coupled, label)
        self.h = h
        self.scaling = scaling
        self.__des = None

    def callback(self):
        """
        Return the C++ callback associated with this Sregion.
        """
        return self.des

    @property
    def h(self):
        """
        Dominance (h)
        """
        return self.__h

    @h.setter
    def h(self, h):
        if math.isinf(h):
            raise ValueError("fwdpy11.Sregion: h not finite")
        if math.isnan(h):
            raise ValueError("fwdpy11.Segion: h not a number")
        self.__h = float(h)

    @property
    def scaling(self):
        """
        Scaling term for the DFE
        """
        return self.__scaling

    @scaling.setter
    def scaling(self, scaling):
        import numpy as np
        self.__scaling = np.uint32(scaling)

    @property
    def des(self):
        """
        Distribution of effect sizes.
        """
        import fwdpy11.fwdpp_extensions
        return getattr(fwdpy11.fwdpp_extensions, self.__des[0])(*self.__des[1])

    @des.setter
    def des(self, des):
        try:
            import fwdpy11.fwdpp_extensions
            if callable(getattr(fwdpy11.fwdpp_extensions, des[0])) is False:
                raise ValueError(
                    "callable required for distribution of effect sizes")
            getattr(fwdpy11.fwdpp_extensions, des[0])(*des[1])
        except:
            raise
        self.__des = des

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

    See :func:`fwdpy11.wright_fisher.evolve` for how this
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
        super(GammaS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)
        self.mean = float(mean)
        self.shape = float(shape)
        self.des = ("makeGammaSH", (self.mean,
                                    self.shape, self.h, self.scaling))

    @property
    def mean(self):
        """
        Mean of the distribution.
        """
        return self.__mean

    @mean.setter
    def mean(self, mean):
        if math.isinf(mean):
            raise ValueError("fwdpy11.GammaS: mean not finite")
        if math.isnan(mean):
            raise ValueError("fwdpy11.GammaS: mean not a number")
        self.__mean = float(mean)

    @property
    def shape(self):
        """
        Shape parameter.
        """
        return self.__shape

    @shape.setter
    def shape(self, shape):
        if math.isinf(shape):
            raise ValueError("fwdpy11.GammaS: shape not finite")
        if math.isnan(shape):
            raise ValueError("fwdpy11.GammaS: shape not a number")
        self.__shape = float(shape)

    def __repr__(self):
        x = 'GammaS(beg=%s,end=%s,weight=%s,'
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

    See :func:`fwdpy11.wright_fisher.evolve` for how this class may be used to
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
        super(ConstantS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)
        self.s = s
        self.des = ("makeConstantSH", (self.s, self.h, self.scaling))

    @property
    def s(self):
        """
        Effect size.
        """
        return self.__s

    @s.setter
    def s(self, s):
        if math.isinf(s):
            raise ValueError("fwdpy11.ConstantS: s not finite")
        if math.isnan(s):
            raise ValueError("fwdpy11.ConstantS: s not a number")
        self.__s = float(s)

    def __repr__(self):
        x = 'ConstantS(beg=%s,end=%s,weight=%s,'
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

    See :func:`fwdpy11.wright_fisher.evolve` for how this
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
            #s is uniform on [-1, 0]
            uniformS = fwdpy11.UniformS(0,1,1,-1,0,0)
        """
        super(UniformS, self).__init__(beg, end, weight, h, coupled, scaling)
        self.lo = lo
        self.hi = hi
        self.des = ("makeUniformSH", (self.lo, self.hi, self.h, self.scaling))

    @property
    def lo(self):
        """
        Lo value for effect size.
        """
        return self.__lo

    @lo.setter
    def lo(self, lo):
        if math.isinf(lo):
            raise ValueError("fwdpy11.UniformS: lo not finite")
        if math.isnan(lo):
            raise ValueError("fwdpy11.UniformS: lo not a number")
        self.__lo = float(lo)

    @property
    def hi(self):
        """
        Hi value for effect size.
        """
        return self.__hi

    @hi.setter
    def hi(self, hi):
        if math.isinf(hi):
            raise ValueError("fwdpy11.UniformS: hi not finite")
        if math.isnan(hi):
            raise ValueError("fwdpy11.UniformS: hi not a number")
        self.__hi = float(hi)

    def __repr__(self):
        x = 'UniformS(beg=%s,end=%s,weight=%s,'
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

    See :func:`fwdpy11.wright_fisher.evolve` for how this class may be used to
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
            expS = fwdpy11.ExpS(0,1,1,-0.1,0)
        """
        super(ExpS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)
        self.mean = mean
        self.des = ("makeExpSH", (self.mean, self.h, self.scaling))

    @property
    def mean(self):
        """
        Mean of the distribution
        """
        return self.__mean

    @mean.setter
    def mean(self, mean):
        if math.isinf(mean):
            raise ValueError("fwdpy11.ExpS: mean not finite")
        if math.isnan(mean):
            raise ValueError("fwdpy11.ExpS: mean not a number")
        self.__mean = float(mean)

    def __repr__(self):
        x = 'ExpS(beg=%s,end=%s,weight=%s,'
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

    See :func:`fwdpy11.wright_fisher.evolve` for how this class may be used to
    parameterize a simulation
    """

    def __init__(self, beg, end, weight, sd, h=1.0, coupled=True, label=0, scaling=1):
        """
        Constructor

        :param beg: the beginning of the region
        :param end: the end of the region
        :param weight: the weight to assign
        :param sd: standard deviation of effect sizes 
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
            gaussianS = fwdpy11.GaussianS(0,1,1,0.1,1)
        """
        super(GaussianS, self).__init__(
            beg, end, weight, h, coupled, label, scaling)
        self.sd = sd
        self.des = ("makeGaussianSH", (self.sd, self.h, self.scaling))

    @property
    def sd(self):
        """
        Standard deviation of the distribution
        """
        return self.__sd

    @sd.setter
    def sd(self, sd):
        if math.isinf(sd):
            raise ValueError("fwdpy11.GaussianS: sd not finite")
        if math.isnan(sd):
            raise ValueError("fwdpy11.GaussianS: sd not a number")
        self.__sd = float(sd)

    def __repr__(self):
        x = 'GaussianS(beg=%s,end=%s,weight=%s,'
        x += 'sd=%s,h=%s,coupled=%s,label=%s,scaling=%s)'
        return x % (self.b, self.e, self.w, self.sd, self.h, self.c, self.l, self.scaling)
