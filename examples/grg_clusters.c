// KMB 2006 Dec 18
// look at mean size of largest cluster of GRG(n,r)
// make grg_clusters && ./grg_clusters | p -x 0 10 -y 0 1

#include <math.h>
#include "vn_graph.h"

int f(graph_t g, int i, void *cd) {
  histogram_t h=*(histogram_t*)cd;
  histogram_add(h,graph_max_cluster(g));
  //return 0; // done
  //fprintf(stdout,"i=%d\n",i);
  if (i<20) return 1; // not done
  if (histogram_mean_done(h,1.0e-2)) return 0; // done
  //if (i>1000) return 0; // done
  return 1; // continue
}

int main(int argc, char* argv[]) {
  unsigned int n=1000,ni=50;
  int i;
  double r,d,dmin=0,dmax=10;
  histogram_t h;
  if (argc>1) n=atoi(argv[1]);    // nodes
  if (argc>2) dmin=atof(argv[2]); // minimum average degree
  if (argc>3) dmax=atof(argv[3]); // maximum average degree
  h=histogram_new("largest cluster of grg");
  for (i=0; i<ni; i++) {
    d=dmin+(dmax-dmin)*i/ni; // average degree
    histogram_empty(h);
    r=sqrt(d/M_PI/(n-1));
    graph_grg_torus_iterator(n,r,f,(void*)&h);
    fprintf(stdout,"%g %g\n",d,histogram_mean(h)/n);
    fflush(stdout);
  }
  histogram_clear(h);
  return 0;
}

