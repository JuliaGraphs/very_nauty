// LLS 2006 Dec 01

#include <math.h>
#include "vn_graph.h"

/*
int graph_largest_cluster(graph_t g) {
  int i,j,k,m=-1,n=graph_nclusters(g);
  for (i=0; i<n; i++) {
    k=0;
    for (j=0; j<nnodes(g); j++) if (cluster(g,j)==i) k++;
    if (k>m) m=k;
  }
  return m;
}

int graph_min_degree(graph_t g) {
  int i,d,k=g->nnodes;
  for (i=0; i<g->nnodes; i++) if ((d=graph_node_degree(g,i))<k) k=d;
  return k;
}
*/

int main(int argc, char* argv[]) {
  int i,j,k,n=1000;
  double r;
  graph_t g;
  histogram_t h;
  if (argc>1) n=atoi(argv[1]);
  g=graph_new(n);
  h=histogram_new("largest cluster");
  for (j=0; j<400; j++) {
    histogram_empty(h);
    r=0.0+0.005*j;
    i=0;
    do {
      graph_grg_torus(g,r);
      k=graph_max_cluster(g);
      //k=graph_min_degree(g);
      //k=graph_max_degree(g);
      //fprintf(stderr,"%d\n",k);
      //fprintf(stdout,"%g %g\n",3.14*a*a*100,(double)k/n);
      //k=graph_nclusters(g);
      histogram_add(h,k);
    } while(++i<50 || !histogram_mean_done(h,1.0e-2));
    //histogram_show(h);
    fprintf(stdout,"%g %g\n",3.14*r*r*10000,histogram_mean(h)/n);
  }
  histogram_clear(h);
  graph_clear(g);
  return 0;
}

