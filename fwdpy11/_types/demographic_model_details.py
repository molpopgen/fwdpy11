import typing

import attr

from ..class_decorators import (attr_add_asblack, attr_class_pickle,
                                attr_class_to_from_dict,
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
    citation: typing.Optional[typing.Union[DemographicModelCitation, typing.Dict]]
    metadata: typing.Optional[object] = None
