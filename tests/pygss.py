import math

import attr

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
