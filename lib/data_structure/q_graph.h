#include "graph_access.h"
#include "definitions.h"
#include <vector>

class q_graph{
private:
  std::vector<Node> nodeList;
  std::vector<Edge> edgeList;
  NodeID num_nodes;
  EdgeID num_edges;
public:
  q_graph(graph_access & G){
    num_nodes = G.number_of_nodes();
    num_edges = G.number_of_edges();
    nodeList = G.get_node_list();
    edgeList = G.get_edge_list();
  }
  
};
