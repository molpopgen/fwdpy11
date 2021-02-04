#include <sstream>
#include <fwdpy11/types/DiploidPopulation.hpp>

static std::string
strip_unix_path(const std::string file)
{
    auto pos = file.find_last_of('/');
    if (pos != std::string::npos)
        {
            return std::string(begin(file) + static_cast<std::ptrdiff_t>(pos) + 1,
                               end(file));
        }
    return file;
}

void
check_mutation_table_consistency_with_count_vectors(const fwdpy11::Population& pop,
                                                    std::string file, int line)
{
    for (auto& mr : pop.tables->mutations)
        {
            if (pop.mcounts[mr.key] + pop.mcounts_from_preserved_nodes[mr.key] == 0)
                {
                    std::ostringstream msg;
                    msg << "mutation table is inconsistent with count vectors: "
                        << strip_unix_path(file) << ", line " << line;
                    throw std::runtime_error(msg.str());
                }
        }
}
