// KMB 2006 Jul 27

#include <math.h>
#include "vn_graph.h"

int graph_largest_cluster(graph_t g) {
  int i,j,k,m=-1,n=graph_nclusters(g);
  for (i=0; i<n; i++) {
    k=0;
    for (j=0; j<nnodes(g); j++) if (cluster(g,j)==i) k++;
    if (k>m) m=k;
  }
  return m;
}

int main(int argc, char* argv[]) {
  int i,j,k,n=500;
  double a,r;
  graph_t g;
  histogram_t h;
  if (argc>1) n=atoi(argv[1]);
  g=graph_new(n);
  for (j=0; j<100; j++) {
    h=histogram_new("largest cluster");
    a=0.0+0.01*j;
    r=a*sqrt(log(n)/n);
    i=0;
    do {
      graph_grg_torus(g,r);
      k=graph_largest_cluster(g);
      //fprintf(stderr,"%d\n",k);
      histogram_add(h,k);
    } while(++i<50 || !histogram_mean_done(h,1.0e-2));
    //histogram_show(h);
    fprintf(stdout,"%g %g\n",a,histogram_mean(h)/n);
    histogram_clear(h);
  }
  graph_clear(g);
  return 0;
}

