// KMB 2006 Nov 22
// process geng output (graph6 or sparse6) and compute edge chromatic number
// e.g.
//   geng 8 | ./vn_graph_chi_dash
//   genrang -P200/1000 20 1000 | ./vn_graph_chi

#include "vn_graph.h"

int chi_dash(graph_t g) {
  return graph_edge_chromatic_number(g,0);
}

int main(int argc, char* argv[]) {
  int i,mx;
  //double x;
  histogram_t h=graph_geng_reader(stdin,chi_dash,"graph_chi_dash");
  histogram_write(stdout,h);
  //x=0.50;
  //fprintf(stderr,"%.2f%% quantile= %g\n",x,histogram_quantile(h,x));
  //fprintf(stderr,"mean=%g\n",histogram_mean(h));
  mx=histogram_max(h);
  for (i=0; i<=mx; i++) {
    printf("%lu ",histogram_bin(h,i));
  }
  printf("\n");
  return 0;
}
