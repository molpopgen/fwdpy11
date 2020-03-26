#ifndef FWDPY11_SERIALIZATION_BACKWARDS_COMPAT_HPP
#define FWDPY11_SERIALIZATION_BACKWARDS_COMPAT_HPP

#include <cstdint>
#include <algorithm>
#include <fwdpp/io/scalar_serialization.hpp>

namespace fwdpy11
{
    namespace serialization
    {
        namespace backwards_compat
        {
            struct deserialize_old_mutation_layout
            // This is the method used for file format version 5,
            // used through version 0.6.2.
            {
                fwdpp::io::scalar_reader reader;
                deserialize_old_mutation_layout() : reader{}
                {
                }
                template <typename streamtype>
                inline fwdpy11::Mutation
                operator()(streamtype &buffer) const
                {
                    fwdpp::uint_t g;
                    bool neutral = false;
                    double pos, s, h;
                    decltype(fwdpy11::Mutation::xtra) xtra;
                    reader(buffer, &g);
                    reader(buffer, &pos);
                    reader(buffer, &s);
                    reader(buffer, &h);
                    reader(buffer, &xtra);
                    std::size_t ns, nh;
                    reader(buffer, &ns);
                    reader(buffer, &nh);
                    std::vector<double> ss, hs;
                    if (ns)
                        {
                            ss.resize(ns);
                            reader(buffer, ss.data(), ns);
                        }
                    if (nh)
                        {
                            hs.resize(ns);
                            reader(buffer, hs.data(), nh);
                        }
                    if (s == 0.0 && std::all_of(begin(ss), end(ss), [](double d) {
                            return d == 0.0;
                        }))
                        {
                            neutral = true;
                        }
                    return fwdpy11::Mutation(neutral, pos, s, h, g, std::move(ss),
                                             std::move(hs), xtra);
                }
            };

            template <typename mcont_t, typename istreamtype>
            void
            read_mutations(istreamtype &in, mcont_t &mutations)
            {
                std::size_t NMUTS;
                fwdpp::io::scalar_reader()(in, &NMUTS);
                deserialize_old_mutation_layout mr;
                for (fwdpp::uint_t i = 0; i < NMUTS; ++i)
                    {
                        mutations.emplace_back(mr(in));
                    }
            }

            template <typename streamtype, typename poptype>
            inline void
            deserialize_population_details(poptype &pop, streamtype &buffer)
            {
                pop.clear();
                fwdpp::io::scalar_reader reader;
                // Step 0: read N
                reader(buffer, &pop.N);
                backwards_compat::read_mutations(buffer, pop.mutations);
                fwdpp::io::read_haploid_genomes(buffer, pop.haploid_genomes);
                fwdpp::io::read_diploids(buffer, pop.diploids);

                // update the mutation counts
                fwdpp::fwdpp_internal::process_haploid_genomes(
                    pop.haploid_genomes, pop.mutations, pop.mcounts);
                fwdpp::io::read_mutations(buffer, pop.fixations);
                if (!pop.fixations.empty())
                    {
                        pop.fixation_times.resize(pop.fixations.size());
                        reader(buffer, &pop.fixation_times[0], pop.fixations.size());
                    }

                // Finally, fill the lookup table:
                for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                    {
                        if (pop.mcounts[i])
                            pop.mut_lookup.emplace(pop.mutations[i].pos,
                                                   static_cast<fwdpp::uint_t>(i));
                    }
            }

        }
    }
}

#endif
