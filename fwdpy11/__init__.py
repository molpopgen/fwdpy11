from .fwdpp_types import *
from .fwdpy11_types import *
from .fwdpy11_sampling import *
from .regions import *

def get_includes():
    """
    Returns absolute path to location of fwdpy11 headers
    """
    import os,fwdpy11
    return os.path.dirname(fwdpy11.__file__)+'/headers'

def get_fwdpp_includes():
    """
    Returns absolute path to location of the fwdpp headers 
    installed along with fwdpy11.
    """
    return get_includes()+'/fwdpp'
