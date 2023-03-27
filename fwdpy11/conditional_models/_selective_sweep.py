#
# Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
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
import typing

import fwdpy11
from fwdpy11.conditional_models import (ConditionalModelOutput,
                                        NewMutationParameters,
                                        SimulationStatus)

from ._track_added_mutation import _track_added_mutation

__disallowed_kwargs = ["until"]


# NOTE: this entire work flow can impose NASTY side-effects:
#  * We are taking lots of references to s.asdict()tuff in the input params.
#  * In general, the input values can be modified by the simulation.
#  * This includes demographic models, genetic value stuff, etc.
#  * Ideally, we'd subjet it all to a copy.deepcopy, and we should test that...
def _selective_sweep(
    rng: fwdpy11.GSLrng,
    pop: fwdpy11._types.DiploidPopulation,
    params: fwdpy11.ModelParams,
    mutation_parameters: NewMutationParameters,
    stopping_condition: typing.Callable[
        [fwdpy11.DiploidPopulation, int, typing.Tuple[float, float, int]],
        SimulationStatus,
    ],
    **kwargs,
) -> ConditionalModelOutput:
    """
    This function is a wrapper around
    :func:`fwdpy11.conditional_models.track_added_mutation`.

    This function requires a `stopping_condition`.
    If `when` is not given as a keyword argument, it is assumed to be 0.
    The keyword argument `until` is not allowed.
    """
    if "when" not in kwargs or kwargs["when"] is None:
        kwargs["when"] = 0

    for d in __disallowed_kwargs:
        if d in kwargs:
            raise ValueError(f"{d} is not a valid kwarg")

    return _track_added_mutation(
        rng,
        pop,
        params,
        mutation_parameters,
        stopping_condition=stopping_condition,
        **kwargs,
    )
