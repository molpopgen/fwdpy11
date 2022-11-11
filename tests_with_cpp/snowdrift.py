import attr

import fwdpy11.class_decorators
import ll_snowdrift


@fwdpy11.class_decorators.attr_class_to_from_dict
@attr.s()
class DiploidSnowdrift(ll_snowdrift._ll_DiploidSnowdrift):
    """
    This is the user-facing Python representation of
    a simple snowdrift model.
    """

    b1 = attr.ib()
    b2 = attr.ib()
    c1 = attr.ib()
    c2 = attr.ib()
    slope = attr.ib()
    p0 = attr.ib()

    def __attrs_post_init__(self):
        # Need to initialize the C++ base class
        super(DiploidSnowdrift, self).__init__(
            self.b1, self.b2, self.c1, self.c2, self.slope, self.p0
        )

    def __getstate__(self):
        # asdict is provided by the 2nd class
        # decorator (reading from bottom to top)
        # We pickle this class as a tuple.
        # Element 0 contains the kwargs + values
        # needed for __init__ and we send along
        # the current phenytpes from the base class
        return (self.asdict(), self.phenotypes)

    def __setstate__(self, t):
        # update the representation of the Python/attrs
        # class
        self.__dict__.update(t[0])
        # Failure to properly init the C++ base class
        # will lead to segmentation faults.
        super(DiploidSnowdrift, self).__init__(**t[0])
        # Reset the base class attribute
        self.phenotypes = t[1]
