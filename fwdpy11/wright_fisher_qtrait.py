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
    def __call__(self,g,e):
        devsq=pow((g+e)-self.O,2)
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
        self.optima = optima.copy()
        self.env = self.optima.pop(0)
    def __call__(self,g,e):
        devsq=pow((g+e)-self.env[1],2)
        return math.exp(-devsq/(2.0*self.env[2]))
    def update(self,generation):
        """
        Update the fitness model conditions.

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
        if(type(rng) is fwdpy11.GSLrng) is False:
            raise ValueError("rng must be a fwdpy11.GSLrng")
        self.sd=sd
        self.mean=mean
        self.rng=rng

    def __call__(self,parent1,parent2):
        return self.mean+fwdpy11.gsl_random.gsl_ran_gaussian_ziggurat(self.rng,self.sd) 

##END LINE NUMBER EMBARGO##

def _evolve_slocus(rng,pop,params,recorder = None):
    import warnings
    #Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        params.validate()
    from .internal import makeMutationRegions,makeRecombinationRegions
    from functools import partial
    mm=makeMutationRegions(params.nregions,params.sregions)
    rm=makeRecombinationRegions(params.recregions)

    noise = None
    if params.noise is None:
        noise = GaussianNoise(rng,0.)
    else:
        noise = params.noise
    updater = None
    noise_updater = None
    if hasattr(params.trait2w,'update'):
        updater = partial(type(params.trait2w).update,params.trait2w)
    if hasattr(noise,'update'):
        noise_updater = partial(type(noise).update,noise)
    if recorder is None:
        from fwdpy11.temporal_samplers import RecordNothing
        recorder = RecordNothing()

    evolve_singlepop_regions_qtrait_cpp(rng,pop,params.demography,
            params.mutrate_n,params.mutrate_s,params.recrate,
            mm,rm,
            params.gvalue,recorder,params.pself,params.trait2w,updater,
            params.noise,noise_updater)

def _evolve_mlocus(rng,pop,params,recorder=None):
    import warnings
    #Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        params.validate()
    
    from .internal import makeMutationRegions,makeRecombinationRegions
    from functools import partial

    noise = None
    if params.noise is None:
        noise = GaussianNoise(rng,0.)
    else:
        noise = params.noise
    mm=[makeMutationRegions(i,j) for i,j in zip(params.nregions,params.sregions)]
    rm=[makeRecombinationRegions(i) for i in params.recregions]
    updater = None
    noise_updater = None
    if hasattr(params.trait2w,'update'):
        updater = partial(type(params.trait2w).update,params.trait2w)
    if hasattr(noise,'update'):
        noise_updater = partial(type(noise).update,noise)
    if recorder is None:
        from fwdpy11.temporal_samplers import RecordNothing
        recorder = RecordNothing()
    evolve_qtrait_mloc_regions_cpp(rng,pop,params.demography,
            params.mutrates_n,params.mutrates_s,params.recrates,
            mm,rm,params.interlocus,
            params.gvalue,
            recorder,
            params.pself,
            params.aggregator,
            params.trait2w,
            updater,noise,noise_updater)

def evolve(rng,pop,params,recorder=None):
    """
    Evolve a quantitative trait.

    This function evolves a class of models where traits are calculated 
    and then mapped to fitness.  See the sections of the manual on 
    simulating quantitative traits.

    :param rng: An instance of :class:`fwdpy11.fwdpy11_types.GSLrng`.
    :param pop: An instance of :class:`fwdpy11.fwdpy11_types.Spop` or :class:`fwdpy11.fwdpy11_types.MlocusPop`.
    :param params: An instance of :class:`fwdpy11.model_params.SlocusParamsQ` or :class:`fwdpy11.model_params.MlocusParamsQ`.
    :param recorder: (None) A callable to record data from the population.

    .. note::
        If recorder is None, :class:`fwdpy11.temporal_samplers.RecordNothing` will be used.

    .. note::
        If params.noise is None, :class:`fwdpy11.wright_fisher_qtrait.GaussianNoise` will be used with mean and 
        standard deviation both set to zero.

    .. note::
        Please be sure to match your population type to your model parameter type.  For example,
        mixing a single-locus population with multi-locus model parameters will trigger an exception.
    """
    import warnings
    #Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        params.validate()
   
    if recorder is None:
        from fwdpy11.temporal_samplers import RecordNothing
        recorder = RecordNothing()
    try:
        _evolve_slocus(rng,pop,params,recorder)
    except:
        _evolve_mlocus(rng,pop,params,recorder)
