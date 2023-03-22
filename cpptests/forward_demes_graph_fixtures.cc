#include "forward_demes_graph_fixtures.hpp"
#include <unistd.h>

const char* single_deme_model = R"(
description: single deme model
time_units: generations
demes:
 - name: A
   epochs:
    - start_size: 100
)";

const char* single_deme_one_size_change = R"(
time_units: generations
demes:
 - name: A
   epochs:
   - start_size: 100
     end_time: 50
   - start_size: 200
)";

const char* two_deme_perpetual_island_model = R"(
time_units: generations
demes:
 - name: A
   epochs:
   - start_size: 100
 - name: B
   epochs:
   - start_size: 100
migrations:
 - demes: [A, B]
   rate: 1e-1
)";

const char* two_deme_perpetual_island_model_with_size_change_and_extinction = R"(
time_units: generations
demes:
 - name: A
   epochs:
   - start_size: 100
     end_time: 25
 - name: B
   epochs:
   - start_size: 100
     end_time: 50
   - start_size: 200
migrations:
 - demes: [A, B]
   rate: 1e-1
)";

const char* two_demes_unequal_merge = R"(
time_units: generations
demes:
 - name: A
   epochs:
    - start_size: 100
      end_time: 25
 - name: B
   epochs:
   - start_size: 75
     end_time: 10
 - name: C
   start_time: 25
   ancestors: [A, B]
   proportions: [0.25, 0.75]
   epochs:
   - start_size: 50
)";

const char* very_recent_pulse_two_generations_ago = R"(
time_units: generations
demes:
 - name: A
   epochs:
   - start_size: 50
 - name: B
   epochs:
   - start_size: 50
pulses:
 - sources: [A]
   dest: B
   proportions: [1.0]
   time: 2
 - sources: [B]
   dest: A
   proportions: [1.0]
   time: 2
)";

const char* extreme_migration_until_one_generation_ago = R"(
time_units: generations
demes:
 - name: A
   epochs:
   - start_size: 1000
 - name: B
   epochs:
   - start_size: 1000
migrations:
 - demes: [A, B]
   rate: 0.5
   end_time: 1
)";

const char* bad_epoch_rounding_02 = R"(
time_units: generations
demes:
- name: bad
  epochs:
  - {end_time: 1.5, start_size: 1}
  - {end_time: 0.4, start_size: 2}
  - {end_time: 0, start_size: 3}
)";

const char* non_integer_start_size = R"(
time_units: generations
demes:
 - name: A
   epochs:
   - start_size: 1000.1
)";

const char* non_integer_end_size = R"(
time_units: generations
demes:
 - name: A
   epochs:
   - start_size: 1000
     end_time: 100
   - end_size: 99.2
)";

const char* linear_size_change = R"(
time_units: generations
demes:
 - name: A
   epochs:
    - start_size: 100
      end_time: 20
    - start_size: 100
      end_size: 200
      end_time: 10
      size_function: linear
    - end_time: 0
      start_size: 55
)";

const char* deme_size_is_one = R"(
time_units: generations
demes:
 - name: A
   epochs:
    - start_size: 1
)";

SingleDemeModel::SingleDemeModel() : yaml(single_deme_model)
{
}

SingleDemeModelOneSizeChange::SingleDemeModelOneSizeChange()
    : yaml(single_deme_one_size_change)
{
}

TwoDemePerpetualIslandModel::TwoDemePerpetualIslandModel()
    : yaml(two_deme_perpetual_island_model)
{
}

TwoDemePerpetualIslandModelWithSizeChangeAndExtinction::
    TwoDemePerpetualIslandModelWithSizeChangeAndExtinction()
    : yaml(two_deme_perpetual_island_model_with_size_change_and_extinction)
{
}

TwoDemesUnequalMerge::TwoDemesUnequalMerge() : yaml(two_demes_unequal_merge)
{
}

VeryRecentPulseTwoGenerationsAgo::VeryRecentPulseTwoGenerationsAgo()
    : yaml(very_recent_pulse_two_generations_ago)
{
}

ExtremeMigrationUntilOneGenerationAgo::ExtremeMigrationUntilOneGenerationAgo()
    : yaml(extreme_migration_until_one_generation_ago)
{
}

BadEpochRounding02::BadEpochRounding02() : yaml(bad_epoch_rounding_02)
{
}

NonIntegerStartSize::NonIntegerStartSize() : yaml(non_integer_start_size)
{
}

NonIntegerEndSize::NonIntegerEndSize() : yaml(non_integer_end_size)
{
}

LinearSizeChange::LinearSizeChange() : yaml(linear_size_change)
{
}

DemeSizeIsOne::DemeSizeIsOne() : yaml(deme_size_is_one)
{
}
