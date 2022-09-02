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

SingleDemeModel::SingleDemeModel() : yaml(single_deme_model)
{
}
