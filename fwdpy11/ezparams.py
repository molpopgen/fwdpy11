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


def mslike(pop, **kwargs):
    """
    Function to establish default parameters
    for a single-locus simulation for standard pop-gen
    modeling scenarios.

    :params pop: An instance of :class:`fwdpy11.DiploidPopulation`
    :params kwargs: Keyword arguments.
    """
    import fwdpy11
    if isinstance(pop, fwdpy11.DiploidPopulation) is False:
        raise ValueError("incorrect pop type: " + str(type(pop)))
    defaults = {'simlen': 10*pop.N,
                'beg': 0.0,
                'end': 1.0,
                'theta': 100.0,
                'pneutral': 1.0,
                'rho': 100.0,
                'dfe': None
                }
    for key, value in kwargs.items():
        if key in defaults:
            defaults[key] = value
    import numpy as np

    params = {'demography': np.array([pop.N]*defaults['simlen'],
              dtype=np.uint32),
              'nregions': [fwdpy11.Region(defaults['beg'],
                           defaults['end'], 1.0)],
              'recregions': [fwdpy11.Region(defaults['beg'],
                             defaults['end'], 1.0)],
              'rates': ((defaults['pneutral']*defaults['theta'])/(4.0*pop.N),
                        ((1.0-defaults['pneutral'])*defaults['theta']) /
                        (4.0*pop.N),
                        defaults['rho']/(4.0*float(pop.N))),
              'gvalue': fwdpy11.Multiplicative(2.0)
              }
    if defaults['dfe'] is None:
        params['sregions'] = []
    else:
        params['sregions'] = [defaults['dfe']]
    return params
