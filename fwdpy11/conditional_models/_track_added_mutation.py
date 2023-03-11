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
import copy
import typing

import attr
import fwdpy11
import numpy as np
from fwdpy11.class_decorators import attr_class_to_from_dict
from fwdpy11.conditional_models import (AddMutationFailure, AlleleCount,
                                        AlleleCountRange, AncientSamplePolicy,
                                        ConditionalModelOutput, EvolveOptions,
                                        FrequencyRange, NewMutationParameters,
                                        OutOfAttempts, SimulationStatus,
                                        _non_negative_value)


@attr_class_to_from_dict
@attr.s(auto_attribs=True)
class _InternalSweepEvolveOptions:
    track_mutation_counts: bool = True


@attr.s(auto_attribs=True)
class _ProgressMonitor:
    index: int = attr.ib(
        validator=[
            attr.validators.instance_of(int),
            _non_negative_value,
        ]
    )
    key: typing.Tuple[float, float, int]
    return_when_stopping_condition_met: bool
    status: SimulationStatus

    def __call__(self, pop: fwdpy11.DiploidPopulation, _) -> bool:
        if (
            self.status.condition_met is True
            and self.return_when_stopping_condition_met is True
        ):
            return True
        return self.status.should_terminate


@attr.s(auto_attribs=True)
class _MutationPresent:
    when: int = attr.ib(
        validator=[attr.validators.instance_of(int), _non_negative_value]
    )
    until: typing.Optional[int] = attr.ib(
        attr.validators.optional(int))  # type: ignore

    def __attrs_post_init__(self):
        if self.until is None:
            raise ValueError(
                "until cannot be None if stopping_condition is None")
        if self.until is not None and self.until <= self.when:
            raise ValueError("until must be > when")

    def is_fixed(self, pop: fwdpy11.DiploidPopulation, key: tuple) -> bool:
        for f in pop.fixations:
            if f.key == key:
                return True
        return False

    def __call__(
        self, pop, index: int, key: tuple
    ) -> SimulationStatus:
        if pop.generation == self.until:
            if pop.mutations[index].key != key and self.is_fixed(pop, key) is False:
                return SimulationStatus(True, False)
            if pop.mcounts[index] > 0:
                return SimulationStatus(False, True)
            if self.is_fixed(pop, key):
                return SimulationStatus(False, True)
        return SimulationStatus(False, False)


@attr.s(auto_attribs=True)
class _Recorder:
    when: typing.Optional[int]
    until: typing.Optional[int]
    criterion: typing.Callable[
        [fwdpy11.DiploidPopulation, int, typing.Tuple[float, float, int]],
        SimulationStatus,
    ]
    monitor: _ProgressMonitor
    sampling_policy: AncientSamplePolicy

    def __attrs_post_init__(self):
        if self.sampling_policy != AncientSamplePolicy.NEVER:
            if self.when is None or self.when < 0:
                raise ValueError(
                    f"sampling policy is {self.sampling_policy} and when is {self.when}"
                )

    def __call__(self, pop, sampler) -> None:
        if self.monitor.status.condition_met is False:
            self.monitor.status = self.criterion(
                pop, self.monitor.index, self.monitor.key
            )
        else:
            return

        x = (
            self.sampling_policy == AncientSamplePolicy.DURATION
            and pop.generation >= self.when
            and (self.until is None or pop.generation <= self.until)
        )

        y = (
            self.sampling_policy == AncientSamplePolicy.COMPLETION
            and self.monitor.status.condition_met is True
            and (self.until is None or pop.generation == self.until)
        )

        if x is True or y is True:
            # Record all alive individuals
            sampler.assign(np.arange(pop.N, dtype=np.uint32))


def _integer_count_details(
    minimum: int,
    maximum: int,
    deme: typing.Optional[int],
    pop: fwdpy11.DiploidPopulation,
):
    if deme is None:
        for c in [minimum, maximum]:
            if c >= 2 * pop.N:
                raise ValueError(
                    f"count {c} "
                    f"is not compatible with total population size of {pop.N}"
                )
    else:
        deme_sizes = pop.deme_sizes(as_dict=True)
        if deme not in deme_sizes:
            raise ValueError(
                f"deme {deme} not present "
                f"at time {pop.generation}. Deme sizes are {deme_sizes}"
            )
        for c in [minimum, maximum]:
            if c >= 2 * pop.N:
                raise ValueError(
                    f"count {c} "
                    f"is not compatible with total population size of {pop.N}"
                )
    return range(minimum, maximum + 1)


def _get_allele_count_range(
    pop: fwdpy11.DiploidPopulation, mutation_parameters: NewMutationParameters
):
    if isinstance(mutation_parameters.frequency, AlleleCount):
        return _integer_count_details(
            mutation_parameters.frequency.count,
            mutation_parameters.frequency.count,
            mutation_parameters.deme,
            pop,
        )
    elif isinstance(mutation_parameters.frequency, AlleleCountRange):
        return _integer_count_details(
            mutation_parameters.frequency.minimum,
            mutation_parameters.frequency.maximum,
            mutation_parameters.deme,
            pop,
        )
    elif isinstance(mutation_parameters.frequency, FrequencyRange):
        if mutation_parameters.deme is None:
            lo = int(np.ceil(mutation_parameters.frequency.minimum * pop.N))
            hi = int(np.floor(mutation_parameters.frequency.minimum * pop.N))
        else:
            deme_sizes = pop.deme_sizes(as_dict=True)
            if mutation_parameters.deme not in deme_sizes:
                raise ValueError(
                    f"deme {mutation_parameters.deme} not present "
                    "at time {pop.generation}"
                )
            lo = int(
                np.ceil(
                    mutation_parameters.frequency.minimum
                    * deme_sizes[mutation_parameters.deme]
                )
            )
            hi = int(
                np.floor(
                    mutation_parameters.frequency.minimum
                    * deme_sizes[mutation_parameters.deme]
                )
            )
        return _integer_count_details(lo, hi, mutation_parameters.deme, pop)
    else:
        raise TypeError(
            f"unsupported type {type(mutation_parameters.frequency)}")


def _copy_pop_and_add_mutation(
    rng,
    pop,
    params,
    mutation_parameters,
    when,
    evolvets_options,
    sampling_policy,
):
    idx = None
    final_count = None
    pcopy = None
    out_params = None
    if when is None or when == 0:
        count_range = _get_allele_count_range(pop, mutation_parameters)

        for c in count_range:
            pcopy = copy.deepcopy(pop)
            _mutation_params = {
                "window": (
                    mutation_parameters.position.left,
                    mutation_parameters.position.right,
                ),
                "ndescendants": c,
                "data": mutation_parameters.data,
                "deme": mutation_parameters.deme,
            }
            idx = pcopy.add_mutation(rng, **_mutation_params)
            if idx is not None:
                final_count = c
                break
        if idx is None:
            raise AddMutationFailure("failed to add mutation")
        out_params = copy.deepcopy(params)
    else:
        if when < 0:
            raise ValueError(f"when must be >= 0, got {when}")

        pcopy = copy.deepcopy(pop)
        pre_sweep_pdict = {k: copy.deepcopy(v)
                           for k, v in params.asdict().items()}
        pre_sweep_pdict["simlen"] = when
        pre_sweep_params = fwdpy11.ModelParams(**pre_sweep_pdict)
        fwdpy11.evolvets(rng, pcopy, pre_sweep_params,
                         **evolvets_options.asdict())

        count_range = _get_allele_count_range(pcopy, mutation_parameters)
        for c in count_range:
            _mutation_params = {
                "window": (
                    mutation_parameters.position.left,
                    mutation_parameters.position.right,
                ),
                "ndescendants": c,
                "data": mutation_parameters.data,
                "deme": mutation_parameters.deme,
            }
            idx = pcopy.add_mutation(rng, **_mutation_params)
            if idx is not None:
                final_count = c
                break

        if idx is None:
            raise AddMutationFailure("failed to add mutation")
        sweep_pdict = {k: v for k, v in pre_sweep_params.asdict().items()}
        sweep_pdict["simlen"] = params.simlen - pcopy.generation
        out_params = fwdpy11.ModelParams(**sweep_pdict)

    return pcopy, idx, final_count, out_params


def _track_added_mutation(
    rng: fwdpy11.GSLrng,
    pop: fwdpy11._types.DiploidPopulation,
    params: fwdpy11.ModelParams,
    mutation_parameters: NewMutationParameters,
    *,
    when: typing.Optional[int] = None,
    until: typing.Optional[int] = None,
    max_attempts: typing.Optional[int] = None,
    sampling_policy: typing.Optional[AncientSamplePolicy] = None,
    stopping_condition=None,
    evolvets_options: typing.Optional[EvolveOptions] = None,
    return_when_stopping_condition_met: bool = False,
) -> ConditionalModelOutput:
    """
    Track the fate of a specific mutation added to the population.

    :param rng: A random number generator
    :type rng: :class:`fwdpy11.GSLrng`
    :param pop: The input population
    :type pop: :class:`fwdpy11.DiploidPopulation`
    :param params: Input simulation parameters
    :type params: :class:`fwdpy11.ModelParams`
    :param mutation_parameters: Details of the new mutation
    :type mutation_parameters: :class:`NewMutationParameters`
    :param when: Time to add the new mutation
    :type when: int
    :param until: Time stop monitoring the new mutation
    :type until: int
    :param max_attempts: Maximum number of attempts to satisfy `stopping_condition`
    :type max_attempts: int
    :param sampling_policy: Enumeration specifying the recording of ancient samples.
    :type sampling_policy: :class:`AncientSamplePolicy`
    :param stopping_condition: An optional callable that specifies a TODO
    :type stopping_condition:  typing.Optional[ typing.Callable[ [fwdpy11.DiploidPopulation, int, typing.Tuple[float, float, int]], SimulationStatus, ]
    :param evolvets_options: Options to :func:`fwdpy11.evolvets`
    :type evolvets_options: :class:`EvolveOptions`
    :param return_when_stopping_condition_met: If `True`, return to calling environment once `stopping_condition` is satisfied.
                                               If `False`, simulate to the end of the model.
    """

    if pop.generation > 0:
        raise ValueError(
            f"the population has been evolved: genration = {pop.generation}"
        )

    if max_attempts is not None:
        if max_attempts < 1:
            raise ValueError(f"max_attempts must be > 0, got {max_attempts}")

    # TODO: all of the stuff below is a hot mess,
    # which would be fixed by the above todo
    # re: arguments to this fn.
    if evolvets_options is None:
        _evolvets_options = EvolveOptions()
    else:
        _evolvets_options = copy.deepcopy(evolvets_options)

    pcopy, idx, ndescendants, local_params = _copy_pop_and_add_mutation(
        rng, pop, params, mutation_parameters, when, _evolvets_options, sampling_policy
    )

    if idx is None:
        raise AddMutationFailure()

    if sampling_policy is None:
        _sampling_policy = AncientSamplePolicy.NEVER
    else:
        _sampling_policy = sampling_policy

    if stopping_condition is None:
        _stopping_condition = _MutationPresent(when, until)  # type: ignore
    else:
        _stopping_condition = stopping_condition

    recorder = _Recorder(
        when,
        until,
        _stopping_condition,
        _ProgressMonitor(
            idx,
            pcopy.mutations[idx].key,
            return_when_stopping_condition_met,
            SimulationStatus(False, False),
        ),
        _sampling_policy,
    )

    internal_options = _InternalSweepEvolveOptions()

    finished = False

    # Rejection-sample our way to glory
    pop_to_return = None

    attempt = 0
    while (
        finished is not True
        and max_attempts is None
        or max_attempts is not None
        and attempt < max_attempts
    ):

        # NOTE: deepcopy and not copy!
        pcopy_loop = copy.deepcopy(pcopy)
        local_params_copy = copy.deepcopy(local_params)
        evolvets_options_copy = copy.deepcopy(_evolvets_options)

        # evolve
        assert evolvets_options_copy is not None
        fwdpy11.evolvets(
            rng,
            pcopy_loop,
            local_params_copy,
            recorder=recorder,
            stopping_criterion=recorder.monitor,
            **evolvets_options_copy.asdict(),  # type: ignore
            **internal_options.asdict(),  # type: ignore
        )

        # The sim ended, so
        # check if condition was satisfied or not
        if recorder.monitor.status.condition_met is True:
            pop_to_return = pcopy_loop
            _evolvets_options = evolvets_options_copy
            local_params = local_params_copy
            finished = True

        attempt += 1

    if attempt == max_attempts:
        raise OutOfAttempts()

    assert pop_to_return is not None
    return ConditionalModelOutput(
        pop=pop_to_return,
        params=local_params,
        mutation_index=idx,
        num_descendant_nodes=ndescendants,
    )
