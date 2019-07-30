#include <fwdpp/ts/marginal_tree_functions/nodes.hpp>

std::vector<fwdpp::ts::TS_NODE_INT>
nodes_preorder(const fwdpp::ts::marginal_tree& m)
{
    return fwdpp::ts::get_nodes(m, fwdpp::ts::nodes_preorder());
}
