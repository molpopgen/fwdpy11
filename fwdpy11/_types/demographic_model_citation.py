import attr

from fwdpy11.class_decorators import (attr_add_asblack, attr_class_pickle,
                                      attr_class_to_from_dict,)


@attr_add_asblack
@attr_class_pickle
@attr_class_to_from_dict
@attr.s(
    kw_only=True, frozen=True, auto_attribs=True, repr_ns="fwdpy11"
)
class DemographicModelCitation(object):
    """
    Citation information for a demographic model

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param DOI: The Digital Object Identifier
    :param full_citation: Something string-like giving the full citation.
    :param metadata: Any additional information that may be needed.


    .. versionadded:: 0.8.0
    """

    DOI: object
    full_citation: object
    metadata: object
