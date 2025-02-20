#include "vn_graph.h"

int main() {
  graph_t g=graph_new(5);
  graph_add_edge(g,0,1);
  graph_add_edge(g,0,2);
  graph_add_edge(g,0,3);
  graph_add_edge(g,0,4);
  graph_add_edge(g,1,4);
  graph_show(g);
  printf("chi=%d\n",graph_chromatic_number(g,0));
  graph_make_dotfile_colored(g,"example_03.dot");
  printf("number of clusters=%d\n",graph_nclusters(g));
  graph_add_node(g);
  printf("nnodes=%d number of clusters=%d\n",nnodes(g),graph_nclusters(g));
  graph_clear(g);
  return 0;
}
