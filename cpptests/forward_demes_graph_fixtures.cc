#include "forward_demes_graph_fixtures.hpp"

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
