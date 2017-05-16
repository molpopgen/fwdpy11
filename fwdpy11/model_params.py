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

"""
    Define classes for parameterizing simulations.

    The classes take either kwargs or an expanded dict as arguments.
    The names of the kwargs are the same as the names of class properties.
    The properties themselves are gettable/settable in the usual way.

    Each class has a validate function that attempts to do a full error-checking
    on the data stored in a class instance.

    The module fwdpy11.ezparams contains some functions for rapidly setting up
    dictionaries of parameters for common cases.

    .. todo::
        The classes should make deep copies of the input so that data 
        from the calling environment are not changed, which can happen 
        in some cases.
"""

import warnings

class ModelParams(object):
    """
    Base class for simulation parameters

    .. versionadded:: 0.1.1
    """
    def __init__(self,**kwargs):
        self.__expected_mutrec_kwargs=['nregions','sregions','recregions']
        self.__expected_demog_kwargs=['demography']
        self.__mutrec_data = {'nregions':[],'sregions':[],'recregions':[]}
        self.__demog_data = {'demography':None}
        self.__fixations = {'prune_selected':None}
        unused_kwargs = []

        for key,value in kwargs.items():
            used = False
            if key in self.__expected_mutrec_kwargs:
                used = True
                self.__mutrec_data[key]=value
            elif key in self.__expected_demog_kwargs:
                used = True
                self.__demog_data[key]=value
            elif key in self.__fixations:
                used = True
                self.__fixations[key]=value
            if used is False:
                unused_kwargs.append(key)

        if len(unused_kwargs) > 0:
            raise ValueError("invalid kwargs: " + ','.join(unused_kwargs))

        if self.__fixations['prune_selected'] is None:
            self.__fixations['prune_selected'] = True

    @property
    def nregions(self):
        return self.__mutrec_data['nregions']

    @property
    def sregions(self):
        return self.__mutrec_data['sregions']

    @property
    def recregions(self):
        return self.__mutrec_data['recregions']

    @property
    def demography(self):
        return self.__demog_data['demography']

    @property 
    def prune_selected(self):
        return self.__fixations['prune_selected']

    @prune_selected.setter
    def prune_selected(self,value):
        self.__fixations['prune_selected'] = bool(value)

    @nregions.setter
    def nregions(self,nregions):
        """
        Set the neutral regions

        :param nregions: The neutral regions.  Type is model-dependant.

        .. note:: See derived class documentation for type details.
        """
        self.__mutrec_data['nregions']=nregions

    @sregions.setter
    def sregions(self,sregions):
        """
        Set the selected regions

        :param sregions: The selected regions.  Type is model-dependant.

        .. note:: See derived class documentation for type details.
        """
        self.__mutrec_data['sregions']=sregions

    @recregions.setter
    def recregions(self,recregions):
        """
        Set the recombination regions.

        :param recregions: The recombination regions.  Type is model-dependant.

        .. note:: See derived class documentation for type details.
        """
        self.__mutrec_data['recregions']=recregions

    @demography.setter
    def demography(self,demog):
        """
        Set the demography

        :param demog: The demographic details, which are model-dependent.

        .. note:: See derived class documentation for type details.
        """
        self.__demog_data['demography']=demog

    def validate(self):
        """
        Error check model params.

        :raises ValueError: Throws ValueError if validation fails.

        .. note:: Demography objects must be validated in derived classes.
        """
        for key,value in self.__mutrec_data.items():
            #We do not check the types 
            #present in each list b/c
            #different modeling scenarios
            #require lists of regions or 
            #lists of lists of regions
            if isinstance(value,list) is False:
                raise ValueError("mutation and recombination data must be lists")

def _validate_single_deme_demography(value):
    import numpy as np
    if any(i<0 for i in value):
        raise ValueError("all population sizes must be > 0")
    if (type(value) is np.ndarray) is False:
        raise ValueError("Type for population size history must be numpy.array")

def _validate_single_locus_rates(data):
    for key,value in data.items():
        if value is None:
            raise ValueError(key + " cannot be None")
        if value < 0.:
            raise ValueError("All mutation and recombination rates must be non-negative.")

class SlocusParams(ModelParams):
    """
    Simulation parameters for single-locus, single-deme simulations.

    .. versionadded:: 0.1.1
    """
    def __set_rates(self,rates):
        try:
            self.__mutrec_data['mutrate_n'] = rates[0]
            self.__mutrec_data['mutrate_s'] = rates[1]
            self.__mutrec_data['recrate'] = rates[2]
            _validate_single_locus_rates(self.__mutrec_data)
        except:
            try:
                if isinstance(rates,dict):
                    for key,value in rates.items():
                        if key in self.__expected_mutrec_kwargs:
                            self.__mutrec_data[key]=value
                        else:
                            raise ValueError("invalid key/value pair for setting rates: ",key,','+str(value))
                else:
                    raise
            except:
                raise
        _validate_single_locus_rates(self.__mutrec_data)

    
    def __init__(self,**kwargs):
        self.__expected_mutrec_kwargs=['mutrate_n','mutrate_s','recrate']
        self.__expected_gvalue_kwargs=['gvalue']
        self.__expected_selfing_kwargs=['pself']
        self.__expected_demog_kwargs=['demography']
        self.__mutrec_data = {};
        self.__gvalue_data = {'gvalue':None}
        self.__selfing_data = {'pself':0.0}

        for i in self.__expected_mutrec_kwargs:
            #Set default mut and rec rates to 0.0
            self.__mutrec_data[i] = 0.0

        new_kwargs = {}
        for key,value in kwargs.items():
            used = False
            if key == 'rates':
                used = True
                try:
                    self.__set_rates(value)
                except:
                    raise
            if key in self.__expected_mutrec_kwargs:
                used = True
                if (type(value) is float) is False:
                    raise ValueError("mutation and recombination rates must be floats")
                self.__mutrec_data[key]=value
            if key in self.__expected_gvalue_kwargs:
                used = True
                self.__gvalue_data[key] = value
            if key in self.__expected_selfing_kwargs:
                used = True
                if (type(value) is float) is False:
                    raise ValueError("selfing probability must be a float")
                self.__selfing_data[key] = value
            if key in self.__expected_demog_kwargs:
                _validate_single_deme_demography(value)
            if used is False:
                new_kwargs[key]=value

        #Set a default genetic value model,
        #which is for "pop-gen".
        #We have to over-ride this potential default
        #condition in SlocusParamsQ
        if self.__gvalue_data['gvalue'] is None:
            from fwdpy11.fitness import SlocusMult
            self.__gvalue_data['gvalue'] = SlocusMult(2.0)

        super(SlocusParams,self).__init__(**new_kwargs)
    
    @property
    def mutrate_n(self):
        """
        Read-only access to neutral mutation rate.
        """
        return self.__mutrec_data['mutrate_n']

    @property
    def mutrate_s(self):
        """
        Read-only access to selected mutation rate.
        """
        return self.__mutrec_data['mutrate_s']

    @property
    def recrate(self):
        """
        Read-only access to recombination rate.
        """
        return self.__mutrec_data['recrate']

    @property
    def pself(self):
        """
        Return selfing probability
        """
        return self.__selfing_data['pself']

    @property
    def rates(self):
        """
        Return all rates as a dict.
        """
        return self.__mutrec_data
    
    @ModelParams.nregions.setter
    def nregions(self,nregions):
        """
        Set the neutral regions

        :param nregions: A list of :class:`fwdpy11.regions.Region.
        """
        ModelParams.nregions.fset(self,nregions)

    @ModelParams.sregions.setter
    def sregions(self,sregions):
        """
        Set the selected regions

        :param sregions:  A list of :class:`fwdpy11.regions.Sregion.
        """
        ModelParams.sregions.fset(self,sregions)

    @ModelParams.recregions.setter
    def recregions(self,recregions):
        """
        Set the recombination regions.

        :param recregions: A list of :class:`fwdpy11.regions.Region.
        """
        ModelParams.recregions.fset(self,recregions)

    @ModelParams.demography.setter
    def demography(self,demog):
        """
        Set the demography

        :param demog: A NumPy array reflecting population sizes over time.
        """
        _validate_single_deme_demography(demog)
        ModelParams.demography.fset(self,demog)

    @rates.setter
    def rates(self,rates):
        """
        Set the neutral mutation, selected mutation, and recombination rates.

        :param rates: A list, tuple, or dict.

        List and tuples must contain the three rates (neutral, selected, recombination).

        Dicts must contain them with the keys 'mutrate_n', 'mutrate_s', and 'recrate'.

        :raises ValueError: Raises exception when bad data are encountered.
        """
        self.__set_rates(rates)    

    @property
    def gvalue(self):
        return self.__gvalue_data['gvalue']

    @gvalue.setter
    def gvalue(self,f):
        from fwdpy11.fitness import SlocusFitness
        if isinstance(f,SlocusFitness) is False:
            raise ValueError("invalid genetic value type: " + str(type(f)))
        self.__gvalue_data['gvalue']=f

    def validate(self):
        """
        Sanity-check parameters.
        """
        super(SlocusParams,self).validate()
        from fwdpy11.fitness import SlocusFitness

        _validate_single_locus_rates(self.__mutrec_data)

        if callable(self.gvalue) is False:
            raise ValueError("genetic value function must be callable")

        #If individual rates are nonzero, then
        #corresponding lists of regions cannot be empty:
        checks = {'mutrate_n':ModelParams.nregions.fget(self),
                'mutrate_s':ModelParams.sregions.fget(self),
                'recrate':ModelParams.recregions.fget(self)}
        for key,value in self.__mutrec_data.items():
            if value > 0.0 and len(checks[key]) == 0:
                raise ValueError(key + ' = ' + str(value) + 
                        ' is associated with empty list of regions.')
            #Issue warning as this may reflect a mistake:
            if value == 0. and len(checks[key]) > 0:
                warnings.warn(key + ' = ' + str(value) + 
                    ' has nonzero region length (' 
                    + str(checks[key]) +
                    ').  Is this intentional?')

        if isinstance(self.__gvalue_data['gvalue'],SlocusFitness) is False:
            raise ValueError("invalid genetic value type: " + type(self.__gvalue_data['gvalue']))

        _validate_single_deme_demography(ModelParams.demography.fget(self))

class SlocusParamsQ(SlocusParams):
    """
    Single locus parameter object for quantitative trait simulations

    .. versionadded:: 0.1.1
    """
    def __init__(self,**kwargs):
        self.__qtrait_params = {'trait_to_fitness': None,'noise': None}
        new_kwargs={}
        for key,value in kwargs.items():
            used = False
            if key in self.__qtrait_params:
                self.__qtrait_params[key]=value
                used = True
            if used is False:
                new_kwargs[key]=value

        #If genetic value fxn not defined,
        #the super will assign mult. fitness
        #here, and so we over-ride that:
        if 'gvalue' not in kwargs:
            from fwdpy11.trait_values import SlocusMultTrait
            new_kwargs['gvalue'] = SlocusMultTrait(2.0)

        super(SlocusParamsQ,self).__init__(**new_kwargs)

        if self.__qtrait_params['trait_to_fitness'] is None:
            from fwdpy11.wright_fisher_qtrait import GSS
            self.__qtrait_params['trait_to_fitness'] = GSS(0,0,1)

        ModelParams.prune_selected.fset(self,False)

    @property
    def trait2w(self):
        return self.__qtrait_params['trait_to_fitness']

    @property
    def noise(self):
        return self.__qtrait_params['noise']

    @trait2w.setter
    def trait2w(self,trait_to_fitness):
        if callable(trait_to_fitness) is False:
            raise ValueError("trait_to_fitness must be a callable.")
        self.__qtrait_params['trait_to_fitness']=trait_to_fitness

    @noise.setter
    def noise(self,noise):
        if callable(noise) is False:
            raise ValueError("noise must be a callable.")
        self.__qtrait_params['noise']=noise

    def validate(self):
        super(SlocusParamsQ,self).validate()

        if callable(self.trait2w) is False:
            raise ValueError("trait to fitness mapper must be a callable")
            
        if self.noise is not None and callable(self.noise) is False:
            raise ValueError("noise function must be callable")

        if self.prune_selected is True:
            raise ValueError("invalid value for prune_selected: " + str(self.prune_selected))

        if self.__qtrait_params['trait_to_fitness'] is None:
            raise ValueError("trait to fitness mapping function cannot be None")

def _validate_multilocus_rates(data):
    for key,value in data.items():
        if value is None:
            raise ValueError(key + " cannot be None")
        if isinstance(value,list) is False:
            raise ValueError(key + " must be a list of float")
        if any(i<0. for i in value):
            raise ValueError("all mutation and recombination rates must be non-negative.")
    #All rate lists must be equal in length
    lengths = set([len(value) for key,value in data.items()])
    if len(lengths) > 1:
        raise ValueError("all lists of mutation and recombintion rates must be equal length.")

class MlocusParams(ModelParams):
    """
    Model parameters for multi-locus simulation

    .. versionadded:: 0.1.1
    """
    def __init__(self,**kwargs):
        self.__expected_demog_kwargs=['demography']
        self.__mutrec_data = {'mutrates_n':[],
                'mutrates_s':[],
                'recrates':[]
                }
        self.__interlocus_rec = {'interlocus':[]}
        self.__gvalue_data = {'gvalue':None}
        self.__selfing_data = {'pself':0.0}
        self.__aggregator = {'agg':None}

        new_kwargs = {}
        for key,value in kwargs.items():
            used = False
            if key == 'rates':
                used = True
                try:
                    self.__set_rates(value)
                except:
                    raise
            if key in self.__mutrec_data:
                used = True
                if any([type(i) is float for i in value]) is False:
                    raise ValueError("mutation and recombination rates must be floats")
                self.__mutrec_data[key]=value
            if key in self.__gvalue_data:
                used = True
                self.gvalue = value
            if key in self.__selfing_data:
                used = True
                if (type(value) is float) is False:
                    raise ValueError("selfing probability must be a float")
                self.__selfing_data[key] = value
            if key in self.__expected_demog_kwargs:
                _validate_single_deme_demography(value)
            if key in self.__interlocus_rec:
                if any(callable(i) for i in value) is False:
                    raise ValueError("interlocus recombination must contain callables")
                self.__interlocus_rec[key]=value
                used = True
            if key in self.__aggregator:
                self.__aggregator[key]=value
                used=True
            if used is False:
                new_kwargs[key]=value

        if self.__aggregator['agg'] is None:
            from fwdpy11.multilocus import AggMultFitness
            self.__aggregator['agg'] = AggMultFitness()

        super(MlocusParams,self).__init__(**new_kwargs)

        ModelParams.prune_selected.fset(self,False)

        #If no genetic value model is defined,
        #we set a default.
        if self.__gvalue_data['gvalue'] is None and len(self.sregions)>0:
            from fwdpy11.multilocus import MultiLocusGeneticValue
            from fwdpy11.fitness import SlocusMult
            self.__gvalue_data['gvalue'] = MultiLocusGeneticValue([SlocusMult(2.0)]*len(self.sregions))

    @property
    def gvalue(self):
        return self.__gvalue_data['gvalue']

    @gvalue.setter
    def gvalue(self,f):
        from fwdpy11.multilocus import MultiLocusGeneticValue
        from fwdpy11.fitness import SlocusFitness
        if isinstance(f,list) is True: 
            if any(isinstance(i,SlocusFitness) for i in f) is False:
                raise ValueError("all elements in list must be single-locus genetic value objects.")
            self.__gvalue_data['gvalue'] = MultiLocusGeneticValue(f)
        elif isinstance(f,MultiLocusGeneticValue) is True:
            self.__gvalue_data['gvalue']=f
        else:
            raise ValueError("invalid genetic value type: " + str(type(f)))


    @property 
    def aggregator(self):
        return self.__aggregator['agg']

    @aggregator.setter
    def aggregator(self,agg):
        if callable(agg) is False:
            raise ValueError("aggregator must be a callable")
        self.__aggregator['agg'] = agg

    @property
    def nregions(self):
        return self.__mutrec_data['nregions']

    @property
    def sregions(self):
        return self.__mutrec_data['sregions']

    @property
    def recregions(self):
        return self.__mutrec_data['recregions']

    @property
    def demography(self):
        return self.__demog_data['demography']

    @property
    def mutrates_n(self):
        """
        Read-only access to neutral mutation rates.
        """
        return self.__mutrec_data['mutrates_n']

    @property
    def mutrates_s(self):
        """
        Read-only access to selected mutation rates.
        """
        return self.__mutrec_data['mutrates_s']

    @property
    def recrates(self):
        """
        Read-only access to recombination rates.
        """
        return self.__mutrec_data['recrates']

    @property
    def pself(self):
        """
        Return selfing probability
        """
        return self.__selfing_data['pself']

    @property
    def rates(self):
        """
        Return all rates as a dict.
        """
        return self.__mutrec_data

    @rates.setter
    def rates(self,rates):
        """
        Set the neutral mutation, selected mutation, and recombination rates.

        :param rates: A list, tuple, or dict.

        List and tuples must contain lists of the three rates (neutral, selected, recombination).

        Dicts must contain them with the keys 'mutrates_n', 'mutrates_s', and 'recrates'.

        :raises ValueError: Raises exception when bad data are encountered.
        """
        self.__set_rates(rates)

    def __set_rates(self,rates):
        try:
            self.__mutrec_data['mutrates_n'] = rates[0]
            self.__mutrec_data['mutrates_s'] = rates[1]
            self.__mutrec_data['recrates'] = rates[2]
            _validate_multilocus_rates(self.__mutrec_data)
        except:
            try:
                if isinstance(rates,dict):
                    for key,value in rates.items():
                        if key in self.__expected_mutrec_kwargs:
                            self.__mutrec_data[key]=value
                        else:
                            raise ValueError("invalid key/value pair for setting rates: ",key,','+str(value))
                else:
                    raise
            except:
                raise
        _validate_multilocus_rates(self.__mutrec_data)

    @ModelParams.nregions.setter
    def nregions(self,nregions):
        """
        Set the neutral regions

        :param nregions: A list of lists of :class:`fwdpy11.regions.Region.
        """
        ModelParams.nregions.fset(self,nregions)

    @ModelParams.sregions.setter
    def sregions(self,sregions):
        """
        Set the selected regions

        :param sregions:  A list of lists of :class:`fwdpy11.regions.Sregion.
        """
        ModelParam.sregions.fset(self,regions)

    @ModelParams.recregions.setter
    def recregions(self,recregions):
        """
        Set the recombination regions.

        :param recregions: A list of lists of :class:`fwdpy11.regions.Region.
        """
        ModelParams.recregions.fset(self,recregions)

    @ModelParams.demography.setter
    def demography(self,demog):
        """
        Set the demography

        :param demog: A NumPy array reflecting population sizes over time.
        """
        _validate_single_deme_demography(demog)
        ModelParams.demography.fset(self,demog)

    @property
    def interlocus(self):
        return self.__interlocus_rec['interlocus']

    @interlocus.setter
    def interlocus(self,interlocus_rec):
        if any(callable(i) for i in interlocus_rec) is False:
            raise ValueError("interlocus_rec must contain callables.")
        self.__interlocus_rec['interlocus'] = interlocus_rec

    def validate(self):
        """
        Sanity-check the class data
        """
        super(MlocusParams,self).validate()

        _validate_multilocus_rates(self.rates)

        if callable(self.aggregator) is False:
            raise ValueError("aggregator must be callable")

        #Lengths of all regions must be the same
        rlens = set((len(self.nregions),len(self.sregions),len(self.recregions)))
        if len(rlens) > 1:
            raise ValueError("length of all region containers must be equal: " +
                    ','.join([str(i) for i in rlens]))

        #Lengths of all lists of rates must be the same
        ratelens = set([len(value) for key,value in self.rates.items()])
        if len(ratelens) > 1:
            raise ValueError("length of per-locus rate containers must be equal: " +
                    ','.join([str(i) for i in rlens]))

        if rlens != ratelens:
            raise ValueError("region and rate containers of different lengths")

        #For k loci, there must be k-1 interlocus-recombination rates
        nloci = rlens.pop()
        if len(self.interlocus)+1 != nloci:
            raise ValueError("invalid parameters: " +
                    str(nloci) + " loci, but only " +
                    str(len(self.interlocus)) + " between-locus recombination functions.")

        if len(self.gvalue) != nloci:
            raise ValueError("Number of loci does not equal number of functions in genetic value object: " +
                    str(nloci) + " != " + str(len(self.aggregator)))

class MlocusParamsQ(MlocusParams):
    """
    Multi-locus model parameters for quantitative trait simulations

    .. versionadded:: 0.1.1
    """
    def __init__(self,**kwargs):
        new_kwargs={}
        self.__qtrait_params={'trait2w':None,
                'noise':None}

        for key,value in kwargs.items():
            used = False
            if key in self.__qtrait_params:
                self.__qtrait_params[key]=value
                used = True
            if used is False:
                new_kwargs[key]=value

        super(MlocusParamsQ,self).__init__(**new_kwargs)

        #If no genetic value model is defined,
        #we set a default.
        if 'gvalue' not in kwargs and len(self.sregions)>0:
            from fwdpy11.multilocus import MultiLocusGeneticValue
            from fwdpy11.trait_values import SlocusMultTrait
            new_kwargs['gvalue'] = MultiLocusGeneticValue([SlocusMultTrait(2.0)]*len(self.sregions))


    @property
    def trait2w(self):
        return self.__qtrait_params['trait2w']

    @trait2w.setter
    def trait2w(self,trait2w):
        if callable(trait2w) is False:
            raise ValueError("trait2w must be a callable.")
        self.__qtrait_params['trait2w']=trait2w

    @property
    def noise(self):
        return self.__qtrait_params['noise']

    @noise.setter
    def noise(self,noise):
        if callable(noise) is False:
            raise ValueError("noise must be a callable.")
        self.__qtrait_params['noise']=noise


    def validate(self):
        """
        Sanity-check class data
        """
        if callable(self.trait2w) is False:
            raise ValueError("trait2w must be a callable.")

        super(MlocusParamsQ,self).validate()
