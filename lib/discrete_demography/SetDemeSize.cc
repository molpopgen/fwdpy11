#include <stdexcept>
#include <fwdpy11/discrete_demography/SetDemeSize.hpp>

namespace fwdpy11
{
    namespace discrete_demography
    {
        SetDemeSize::SetDemeSize(std::uint32_t w, std::int32_t d, std::uint32_t n,
                                 bool reset)
            : when(w), deme(d), new_size(n), resets_growth_rate(reset)
        {
            if (deme < 0)
                {
                    throw std::invalid_argument("SetDemeSize:"
                                                " deme must be non-negative");
                }
        }
    }
}
