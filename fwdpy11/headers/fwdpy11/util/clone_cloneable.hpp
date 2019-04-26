#ifndef FWDPY11_UTIL_CLONE_CLONEABLE_HPP
#define FWDPY11_UTIL_CLONE_CLONEABLE_HPP
#include <memory>

namespace fwdpy11
{
    namespace util
    {
        template <typename T>
        inline std::unique_ptr<T>
        clone_cloneable(const std::unique_ptr<T>& g)
        // Much of the API is based on ABC class hierarchies
        // defined in C++.  Python types hold unique_ptrs
        // to these bases on the C++ side.  We often 
        // need to clone the object when doing thing
        // like pickling.  This function abstracts out
        // that op.  Further, note that a 
        // unique_ptr holding nullptr registers as 
        // None in python, which is a very nice feature.
        {
            if (g == nullptr)
                {
                    return std::unique_ptr<T>(nullptr);
                }
            return g->clone();
        }
    } // namespace util
} // namespace fwdpy11
#endif
