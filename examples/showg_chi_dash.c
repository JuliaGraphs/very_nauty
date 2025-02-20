// KMB 2007 Jan 08
// ./random_planar.py && ./showg_chi_dash < foo.showg

#include "vn_graph.h"

int f(graph_t g) {
  int timeout=10,chi_dash=graph_edge_chromatic_number(g,timeout);
  //fprintf(stderr,"%d ",chi_dash);
  return chi_dash;
}

int main(int argc, char* argv[]) {
  histogram_t h=graph_showg_reader(stdin,f,"chi_dash");
  histogram_show(h);
  fprintf(stderr,"\n");
  return 0;
}
