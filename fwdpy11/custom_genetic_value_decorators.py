def _make_cloner(cls, base, *args):
    def clone(self):
        # create a new object without initializing it
        cloned = cls.__new__(cls)
        # clone C++ state
        base.__init__(cloned, *args)
        # clone Python state
        cloned.__dict__.update(self.__dict__)
        return cloned

    cls.clone = clone
    return cls


class genetic_value_is_trait_default_clone(object):
    def __init__(self, ndim=1):
        self.ndim = ndim

    def __call__(self, cls):
        from fwdpy11 import GeneticValueIsTrait

        return _make_cloner(cls, GeneticValueIsTrait, self.ndim)


def genetic_value_noise_default_clone(cls):
    from fwdpy11 import GeneticValueNoise

    return _make_cloner(cls, GeneticValueNoise)
    return cls


def default_update(cls):
    def update(self, pop):
        pass

    cls.update = update
    return cls
