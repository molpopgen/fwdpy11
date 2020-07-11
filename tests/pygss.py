import math

import attr
import numpy as np

import fwdpy11
import fwdpy11.custom_genetic_value_decorators


@fwdpy11.custom_genetic_value_decorators.default_update
@fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone()
@attr.s()
class PyGSS(fwdpy11.GeneticValueIsTrait):
    opt = attr.ib()
    VS = attr.ib()

    def __attrs_post_init__(self):
        fwdpy11.GeneticValueIsTrait.__init__(self)

    def __call__(self, metadata, genetic_values):
        return math.e ** (-((metadata.g - self.opt) ** 2) / (2 * self.VS))


@fwdpy11.custom_genetic_value_decorators.genetic_value_is_trait_default_clone()
@attr.s()
class PyGSSRandomOptimum(fwdpy11.GeneticValueIsTrait):
    opt = attr.ib()
    VS = attr.ib()

    def __attrs_post_init__(self):
        fwdpy11.GeneticValueIsTrait.__init__(self)
        self.optima = []

    def __call__(self, metadata, genetic_values):
        return math.e ** (-((metadata.g - self.opt) ** 2) / (2 * self.VS))

    def update(self, pop):
        self.opt = np.random.normal(0.0, 1.0, 1,)[0]
        self.optima.append((pop.generation, self.opt))
