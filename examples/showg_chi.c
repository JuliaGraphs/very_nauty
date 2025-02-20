// KMB 2007 Jan 08
// ./random_planar.py 15 100 && ./showg_chi < foo.showg; ./showg_chi_dash < foo.showg
// ./random_planar.py && ./showg_chi < foo.showg

#include "vn_graph.h"

int f(graph_t g) {
  int timeout=10;
  return graph_chromatic_number(g,timeout);
}

int main(int argc, char* argv[]) {
  histogram_t h=graph_showg_reader(stdin,f,"chi");
  histogram_show(h);
  return 0;
}
