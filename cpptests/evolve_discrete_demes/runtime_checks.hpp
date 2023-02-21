#include <string>
#include <fwdpy11/types/Population.hpp>

std::string strip_unix_path(const std::string file);
void check_mutation_table_consistency_with_count_vectors(const fwdpy11::Population& pop,
                                                         std::string file, int line);

