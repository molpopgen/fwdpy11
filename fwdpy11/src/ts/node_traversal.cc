#include <fwdpp/ts/marginal_tree_functions/nodes.hpp>

std::vector<fwdpp::ts::table_index_t>
nodes_preorder(const fwdpp::ts::marginal_tree& m)
{
    return fwdpp::ts::get_nodes(m, fwdpp::ts::nodes_preorder());
}
