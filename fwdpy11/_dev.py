#
# Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
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


def get_includes():
    """
    Returns absolute path to location of fwdpy11 headers
    """
    import os
    import fwdpy11
    return os.path.dirname(fwdpy11.__file__)+'/headers'


def get_fwdpp_includes():
    """
    Returns absolute path to location of the fwdpp headers
    installed along with fwdpy11.
    """
    return get_includes()+'/fwdpp'


def minimal_mako():
    """
    Returns a minimal mako header for compiling
    plugins using cppimport

    .. versionadded:: 0.1.1
    """
    rv = r"""<%
setup_pybind11(cfg)
import fwdpy11 as fp11
cfg['include_dirs'] = [fp11.get_includes(), fp11.get_fwdpp_includes()]
%>"""
    return rv
