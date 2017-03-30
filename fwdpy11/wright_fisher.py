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
from .wfevolve import evolve_singlepop_regions_cpp

def quick_sim(ngens = None):
    """
    A convenience function for rapidly creating a
    :class:`fwdpy11.fwdpy11_types.Spop`

    .. testcode:: 

        import fwdpy11.wright_fisher
        #This will simulate N=1e3 for 10N generations
        pop = fwdpy11.wright_fisher.quick_sim()
        print(pop.N)
        print(pop.generation)

    The output is:

    .. testoutput::

      1000
      10000

    .. note::
        Implemented via a call to :func:`fwdpy11.wright_fisher.evolve`
    """
    from .fwdpy11_types import GSLrng,Spop
    rng = GSLrng(42)
    pop=Spop(1000)
    if ngens is None:
        evolve(rng,pop)
    else:
        import numpy as np
        nlist = np.array([pop.N]*ngens,dtype=np.uint32)
        evolve(rng,pop,nlist)
    return pop

def evolve(rng,pop,popsizes = None,mu_neutral=None,
        mu_selected = None,recrate=None,sregions=None):
    """
    
    .. testcode::

        import fwdpy11 as fp11
        import fwdpy11.wright_fisher as wf
        import numpy as np
        rng = fp11.GSLrng(42)
        p = fp11.Spop(1000)
        wf.evolve(rng,p)
        print(p.generation)
        nlist=np.array([5000]*323,dtype=np.uint32)
        wf.evolve(rng,p,nlist)
        print(p.N)
        print(p.generation)
            
    .. testoutput::

      10000
      5000
      10323

    """
    if popsizes is None:
        import numpy as np
        popsizes = np.array([pop.N]*10*pop.N,dtype=np.uint32)
    if mu_neutral is None:
        mu_neutral = 100./float(4*pop.N)
    if mu_selected is None:
        mu_selected = 10./float(4*pop.N)
    if recrate is None:
        recrate = 100./float(4*pop.N)
    if sregions is None:
        from .regions import ExpS
        sregions = [ExpS(0,1,1,-0.1,1.0)]
    from .regions import Region
    nr=[Region(0,1,1)]
    return evolve_regions(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,nr,sregions,
            nr)

def evolve_regions(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,
        selfing_rate = 0.):
    """
    Evolve a single deme according to a Wright-Fisher life cycle 
    with arbitrary changes in population size and a temporal sampler.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param pop: A :class:`fwdpy11.fwdpy11_types.Spop`
    :param popsizes: A 1d NumPy array representing population sizes over time.
    :param mu_neutral: The neutral mutation rate (per gamete, per generation)
    :param mu_selected: The selected mutation rate (per gamete, per generation)
    :param recrate: The recombination reate (per diploid, per generation)
    :param nregions: A list of :class:`fwdpy11.regions.Region`.
    :param sregions: A list of :class:`fwdpy11.regions.Sregion`.
    :param recregions: A list of :class:`fwdpy11.regions.Region`.
    :param recorder: A callable to record data from the population.
    :param selfing_rate: (default 0.0) The probability than an individual selfs.

    .. note:: 
        The fitness model will be :class:`fwdpy11.fitness.SpopAdditive` constructed
        with a scaling of 2.0. This function calls 
        :func:`fwdpy11.wright_fisher.evolve_regions_sampler`, passing in a
        :class:`fwdpy11.temporal_samplers.RecordNothing` object.
    """
    from .temporal_samplers import RecordNothing
    recorder=RecordNothing()
    return evolve_regions_sampler(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,nregions,sregions,recregions,
            recorder,selfing_rate)

def evolve_regions_sampler(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,
        recorder,selfing_rate = 0.):
    """
    Evolve a single deme according to a Wright-Fisher life cycle 
    with arbitrary changes in population size and a temporal sampler.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param pop: A :class:`fwdpy11.fwdpy11_types.Spop`
    :param popsizes: A 1d NumPy array representing population sizes over time.
    :param mu_neutral: The neutral mutation rate (per gamete, per generation)
    :param mu_selected: The selected mutation rate (per gamete, per generation)
    :param recrate: The recombination reate (per diploid, per generation)
    :param nregions: A list of :class:`fwdpy11.regions.Region`.
    :param sregions: A list of :class:`fwdpy11.regions.Sregion`.
    :param recregions: A list of :class:`fwdpy11.regions.Region`.
    :param recorder: A callable to record data from the population.
    :param selfing_rate: (default 0.0) The probability than an individual selfs.

    .. note:: 
        The fitness model will be :class:`fwdpy11.fitness.SpopAdditive` constructed
        with a scaling of 2.0.
    """
    from .fitness import SpopAdditive
    fitness = SpopAdditive(2.0)
    return evolve_regions_sampler_fitness(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,nregions,sregions,recregions,fitness,
            recorder,selfing_rate)

def evolve_regions_sampler_fitness(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,fitness,
        recorder,selfing_rate = 0.):
    """
    Evolve a single deme according to a Wright-Fisher life cycle 
    with arbitrary changes in population size, a specified fitness model,
    and a temporal sampler.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param pop: A :class:`fwdpy11.fwdpy11_types.Spop`
    :param popsizes: A 1d NumPy array representing population sizes over time.
    :param mu_neutral: The neutral mutation rate (per gamete, per generation)
    :param mu_selected: The selected mutation rate (per gamete, per generation)
    :param recrate: The recombination reate (per diploid, per generation)
    :param nregions: A list of :class:`fwdpy11.regions.Region`.
    :param sregions: A list of :class:`fwdpy11.regions.Sregion`.
    :param recregions: A list of :class:`fwdpy11.regions.Region`.
    :param fitness: A :class:`fwdpy11.fitness.SpopFitness`.
    :param recorder: A callable to record data from the population.
    :param selfing_rate: (default 0.0) The probability than an individual selfs.
    """
    from .internal import makeMutationRegions,makeRecombinationRegions
    mm=makeMutationRegions(nregions,sregions)
    rm=makeRecombinationRegions(recregions)
    evolve_singlepop_regions_cpp(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,mm,rm,fitness,recorder,selfing_rate)
