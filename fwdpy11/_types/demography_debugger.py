import fwdpy11
import numpy as np


class DemographyDebugger(object):
    """
    Debug demographic events efficiently.
    """

    def __init__(self, pop, events):
        """
        :param events: The demographic events
        :type events: :class:`fwdpy11.DiscreteDemography`
        """
        self.M = None
        if events.migmatrix is not None:
            self.M = np.copy(events.migmatrix.M)

        self._validate_migration_rate_change_lengths(events)

        self.maxdemes = self._get_maxdemes(pop, events)

        if self.maxdemes < 1:
            raise ValueError("Invalid number of "
                             "demes in simulation: {}".format(self.maxdemes))

        self.current_deme_sizes = np.zeros(self.maxdemes, dtype=np.uint32)
        for i, j in zip(*pop.deme_sizes()):
            self.current_deme_sizes[i] = j

        self.growth_rates = np.array([fwdpy11.NOGROWTH]*self.maxdemes)
        self.growth_onset_times = np.zeros(self.maxdemes, dtype=np.uint32)

    def _get_maxdemes(self, pop, events):
        """
        The maximum number of demes the sim can ever see.
        """
        md = np.array(pop.diploid_metadata, copy=False)
        max_from_md = md['deme'].max()
        max_from_events = -1

        def update_max_from_events(m, e):
            for i in e:
                try:
                    m = max(m, i.deme)
                except AttributeError:
                    m = max(m, i.source)
            return m

        for i in [events.mass_migrations, events.set_growth_rates,
                  events.set_deme_sizes, events.set_selfing_rates,
                  events.set_migration_rates]:
            max_from_events = update_max_from_events(max_from_events, i)

        current_max = max(max_from_md, max_from_events) + 1
        if self.M is None:
            return current_max

        if self.M.shape[0] < current_max:
            raise ValueError("The MigrationMatrix shape, {}, "
                             "does not match the max number of "
                             "demes present in the "
                             "simulation, ".format(self.M.shape, current_max))

        return max(current_max, self.M.shape[0])

    def _validate_migration_rate_change_lengths(self, events):
        """
        Various checks on the lengths of new migration rates.
        """
        if self.M is None and len(events.set_migration_rates) > 0:
            raise ValueError("migration rate changes are "
                             "set but there is no migration matrix")
        for i in events.set_migration_rates:
            if i.deme >= 0:
                if len(i.migrates) != self.M.shape[0]:
                    raise ValueError("Migration rates mismatch")
            else:  # Are replacing the entire matrix
                if len(i.migrates) != len(self.M.flatten()):
                    raise ValueError("Migration rates mismatch")
