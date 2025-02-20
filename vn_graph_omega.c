// KMB 2006 Jul 18
// process geng output (graph6 or sparse6) and compute clique number
// e.g.
//   geng 8 | ./graph_omega
//   genrang -P200/1000 20 1000 | ./graph_omega

#include "vn_graph.h"

inline int chi(graph_t g) {
  return graph_clique_number(g);
}

int main() {
  double x;
  histogram_t h=graph_geng_reader(stdin,chi,"graph_omega");
  histogram_show(h);
  x=0.50;
  printf("%.2f%% quantile= %g\n",x,histogram_quantile(h,x));
  //printf("mean=%g\n",histogram_mean(h));
  return 0;
}
