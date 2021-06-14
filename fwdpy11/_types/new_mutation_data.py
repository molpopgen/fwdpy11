import typing

import attr

from .._fwdpy11 import ll_NewMutationData


@attr.s(auto_attribs=True, frozen=True, repr_ns="fwdpy11", kw_only=True)
class NewMutationData(ll_NewMutationData):
    """
    Data object for :func:`fwdpy11.DiploidPopulation.add_mutation`

    .. versionadded:: 0.16.0

    This class has the following attributes, whose names
    are also `kwargs` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param effect_size: The fixed/univariate effect size.
                        Will fill in :attr:`fwdpy11.Mutation.s`
    :type effect_size: float
    :param dominance: Heterozygous effect of the mutation.
                      Will fill in :attr:`fwdpy11.Mutation.h`
    :type dominance: float
    :param esizes: Data to fill :attr:`fwdpy11.Mutation.esizes`.
                   Default is `None`.
    :type esizes: list[float]
    :param heffects: Data to fill :attr:`fwdpy11.Mutation.heffects`.
                     Default is `None`.
    :type heffects: list[float]
    :param label: Data for :attr:`fwdpy11.Mutation.label`.
                  Default is 0.
    :type label: int
    """

    effect_size: float
    dominance: float
    esizes: typing.Optional[typing.List[float]] = None
    heffects: typing.Optional[typing.List[float]] = None
    label: int = 0

    def __attrs_post_init__(self):
        if self.esizes is None and self.heffects is None:
            super(NewMutationData, self).__init__(
                self.effect_size, self.dominance, [], [], self.label
            )
        else:
            super(NewMutationData, self).__init__(
                self.effect_size, self.dominance, self.esizes, self.heffects, self.label
            )
