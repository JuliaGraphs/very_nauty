// KMB 2006 Jul 18
// process geng output (graph6 or sparse6) and compute chromatic number
// e.g.
//   geng 8 | ./vn_graph_chi
//   genrang -P200/1000 20 1000 | ./vn_graph_chi

#include "vn_graph.h"

int chi(graph_t g) {
  return graph_chromatic_number(g,0);
}

int main() {
  double x;
  histogram_t h=graph_geng_reader(stdin,chi,"graph_chi");
  histogram_show(h);
  x=0.50;
  printf("%.2f%% quantile= %g\n",x,histogram_quantile(h,x));
  printf("  histogram mean=%g\n",histogram_mean(h));
  return 0;
}
