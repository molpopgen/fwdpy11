class genetic_value_is_trait_default_clone(object):
    def __init__(self, ndim=1):
        self.ndim = ndim

    def __call__(self, cls):
        from fwdpy11 import GeneticValueIsTrait

        def clone(me):
            # create a new object without initializing it
            cloned = cls.__new__(cls)
            # clone C++ state
            GeneticValueIsTrait.__init__(cloned, self.ndim)
            # clone Python state
            cloned.__dict__.update(me.__dict__)
            return cloned

        cls.clone = clone
        return cls


def genetic_value_noise_default_clone(cls):
    from fwdpy11 import GeneticValueNoise

    def clone(self):

        # create a new object without initializing it
        cloned = cls.__new__(cls)
        # clone C++ state
        GeneticValueNoise.__init__(cloned)
        # clone Python state
        cloned.__dict__.update(self.__dict__)
        return cloned

    cls.clone = clone
    return cls


def default_update(cls):
    def update(self, pop):
        pass

    cls.update = update
    return cls
