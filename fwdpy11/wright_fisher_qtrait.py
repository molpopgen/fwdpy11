from .wfevolve_qtrait import evolve_singlepop_regions_qtrait_cpp,evolve_qtrait_mloc_regions_cpp
import fwdpy11.gsl_random
import math

#Note: making changes until "END LINE NUMBER EMBARGO"
#is seen below requires changes to the documentation as 
#well!  The fwdpy11 manual includes specific line ranges
#from this file in order to provide examples of trait ->
#fitness and noise functions.  The relevant manual file
#is doc/examples/qtraits.rst

class GSS:
    """
    Gaussian stabilizing selection with a constant optimum.

    This is a callable object.  The __call__ function
    takes a trait value and return a fitness.

    Its parameters are given below.

    :param trait_value: A trait (phenotype) value

    :return: :math:`w=e^{-\\frac{(G-E)^2}{2VS}}`

    :rtype: float
    """
    def __init__(self,VS,O):
        """
        :param VS: 1/VS is intensity of selection against phenotypic deviations from the mean/optimum.
        :param O: The optimum trait value.
        """
        if VS <= 0.:
            raise ValueError("VS > 0 required")
        self.VS=VS
        self.O=O
    def __call__(self,trait_value):
        devsq=pow(trait_value-self.O,2)
        return math.exp(-devsq/(2.0*self.VS))

class GSSmo:
    """
    Gaussian stabilizing selection with moving optimum.

    .. note::
        This class's constructor enforces nothing about the initial 
        values for the model parameters. The reason is for flexibility,
        in that you may run a simulation, and then evolve more later 
        one.  However, it makes sense that, when you start a simulation,
        the first tuple for the constructor refers to the population's
        current generation.
    """
    def __init__(self,optima):
        """
        :param optima: A list of tuples.  Each tuple is (generation,optimum,VS)
        """
        if len(optima) == 0:
            raise ValueError("empty list of optima")
        for oi in optima:
            if isinstance(oi,tuple) is False: 
                raise ValueError("optima must cointain tuples")
            if len(oi) != 3:
                raise ValueError("incorrect tuple length")
            if oi[0] < 0:
                raise ValueError("negative generation not allowed")
            if oi[2] <= 0.0:
                raise ValueError("VS > 0 required")
        self.optima = optima
        self.env = self.optima.pop(0)
    def __call__(self,P):
        devsq=pow(P-self.env[1],2)
        return math.exp(-devsq/(2.0*self.env[2]))
    def update(self,generation):
        """
        Update the fitness model conditins.

        :param generation: the generation in the simulation

        .. note:: this function is called from within the simulation
        """
        if len(self.optima)==0:
            return
        if generation >= self.optima[0][0]:
            self.env = self.optima.pop(0)

class GaussianNoise:
    """
    Gaussian noise for trait values.
    
    Adds :math:`N(\\mu,\\sigma)` to trait values.
    """
    def __init__(self,rng,sd,mean=0.0):
        """
        :param rng: A :class:`fwdpy11.GSLrng`
        :param sd: :math:`\\sigma`
        :param mean: (0.0) :math:`\\mu`
        """
        self.sd=sd
        self.mean=mean
        self.rng=rng

    def __call__(self,g,parent1,parent2):
        return self.mean+fwdpy11.gsl_random.gsl_ran_gaussian_ziggurat(self.rng,self.sd) 

##END LINE NUMBER EMBARGO##

def evolve_regions_sampler_fitness(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,trait_model,
        trait_to_fitness=GSS(1,0.0),recorder=None,noise=None,selfing_rate = 0.):

    """
    Evolve a single deme according to a Wright-Fisher life cycle 
    with arbitrary changes in population size, a specified fitness model,
    and a temporal sampler.

    This function evolves a class of models where traits are calculated 
    and then mapped to fitness.  See the sections of the manual on 
    simulating quantitative traits.

    :param rng: A :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param pop: A :class:`fwdpy11.fwdpy11_types.Spop`
    :param popsizes: A 1d NumPy array representing population sizes over time.
    :param mu_neutral: The neutral mutation rate (per gamete, per generation)
    :param mu_selected: The selected mutation rate (per gamete, per generation)
    :param recrate: The recombination reate (per diploid, per generation)
    :param nregions: A list of :class:`fwdpy11.regions.Region`.
    :param sregions: A list of :class:`fwdpy11.regions.Sregion`.
    :param recregions: A list of :class:`fwdpy11.regions.Region`.
    :param trait_model: A :class:`fwdpy11.fitness.SpopFitness`.
    :param trait_to_fitness: (fwdpy11.wright_fisher_qtrait.GSS(0,0,0)) A callable relating trait value to fitness.
    :param recorder: (None) A callable to record data from the population.
    :param noise: (None) A callable for adding random effects to trait values.


    .. note::
        If recorder is None, :class:`fwdpy11.temporal_samplers.RecordNothing` will be used.

    .. note::
        If noise is None, :class:`fwdpy11.wright_fisher_qtrait.GaussianNoise` will be used with mean and 
        standard deviation both set to zero.
    """
    from .internal import makeMutationRegions,makeRecombinationRegions
    from functools import partial
    if noise is None:
        noise = GaussianNoise(rng,0.)
    mm=makeMutationRegions(nregions,sregions)
    rm=makeRecombinationRegions(recregions)
    updater = None
    noise_updater = None
    if hasattr(trait_to_fitness,'update'):
        updater = partial(type(trait_to_fitness).update,trait_to_fitness)
    if hasattr(noise,'update'):
        noise_updater = partial(type(noise).update,noise)
    if recorder is None:
        from fwdpy11.temporal_samplers import RecordNothing
        recorder = RecordNothing()
    evolve_singlepop_regions_qtrait_cpp(rng,pop,popsizes,mu_neutral,
            mu_selected,recrate,mm,rm,trait_model,recorder,selfing_rate,
            trait_to_fitness,updater,noise,noise_updater)

def evolve_mlocus_regions_sampler_fitness(rng,pop,popsizes,mu_neutral,
        mu_selected,recrate,nregions,sregions,recregions,interlocus_rec,genetic_value_model,
        multilocus_trait_model,
        trait_to_fitness=GSS(1,0.0),recorder=None,noise=None,selfing_rate = 0.):
    from .internal import makeMutationRegions,makeRecombinationRegions
    from functools import partial

    if noise is None:
        noise = GaussianNoise(rng,0.)
    mm=[makeMutationRegions(i,j) for i,j in zip(nregions,sregions)]
    rm=[makeRecombinationRegions(i) for i in recregions]
    updater = None
    noise_updater = None
    if hasattr(trait_to_fitness,'update'):
        updater = partial(type(trait_to_fitness).update,trait_to_fitness)
    if hasattr(noise,'update'):
        noise_updater = partial(type(noise).update,noise)
    if recorder is None:
        from fwdpy11.temporal_samplers import RecordNothing
        recorder = RecordNothing()
    evolve_qtrait_mloc_regions_cpp(rng,pop,popsizes,
            mu_neutral,mu_selected,recrate,
            mm,rm,interlocus_rec,
            genetic_value_model,
            recorder,
            selfing_rate,
            multilocus_trait_model,
            trait_to_fitness,
            updater,noise,noise_updater)
    
