import typing

import attr

from ..class_decorators import (attr_add_asblack, attr_class_pickle,
                                attr_class_to_from_dict_no_recurse)


from .demographic_model_citation import DemographicModelCitation
from .forward_demes_graph import ForwardDemesGraph


@attr_add_asblack
@attr_class_pickle
@attr_class_to_from_dict_no_recurse
@attr.s(
    kw_only=True, frozen=True, auto_attribs=True, repr_ns="fwdpy11"
)
class DemographicModelDetails(object):
    """
    Stores rich information about a demographic model.
    Instances of this class get returned by functions
    generating pre-calculated models.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param model: The demographic model parameters
    :type model: fwdpy11.ForwardDemesGraph
    :param name: The name of the model
    :type name: str
    :param source: The source of the model
    :type source: dict
    :param parameters: The parameters used to generate ``model``
    :type parameters: object
    :param citation: Citation information
    :type citation: dict or fwdpy11.DemographicModelCitation or None
    :param metadata: Optional field for additional info
    :type metadata: object

    .. versionadded:: 0.8.0
    """
    model: ForwardDemesGraph
    name: str
    source: typing.Dict
    parameters: object
    citation: typing.Optional[typing.Union[DemographicModelCitation,
                                           typing.Dict]]
    metadata: typing.Optional[object] = None

    @property
    def initial_sizes_list(self) -> typing.List[int]:
        """
        A list of the nonzero (parental) deme sizes at time 0 of the model
        (thinking forwards in time).

        This property is useful for setting up a
        :class:`fwdpy11.DiploidPopulation`.
        """
        return self.model.initial_sizes

    @property
    def total_simulation_length(self) -> int:
        """
        The total number of generations in the model.

        This value includes any burn-in time and starts
        from time 0 (forwards in time).
        """
        return self.metadata["total_simulation_length"]

    @property
    def deme_labels(self) -> typing.Dict[int, str]:
        """
        A dictionary mapping integer deme IDs to deme names
        (strings).
        """
        return self.metadata["deme_labels"]
