// KMB 2006 Dec 04, 06, 18
// make gnm_clusters
// cmin=1; cmax=6; time nice ./gnm_clusters 10000 $cmin $cmax 1000 | colex 1 3 | p -m0 -x $cmin $cmax -y 0 1
// make gnm_clusters && time ./gnm_clusters 1000000 1 2 100 | colex 3 4 | dm "x1/x2" | histogram | bars | p

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vn_graph.h"
#include "LambertW.c"

double expected(int n, double c) {
  /* Bollob\'as Modern graph theory p 241.
     Typical size of largest component of G(n,m=c*n/2).  */
  double alpha=c-1-log(c);
  if (c<1)  return (log(n)-2.5*log(log(n)))/alpha/n;
  if (c==1) return pow(n,2.0/3)/n;
  return 1.0+LambertW(-c/exp(c))/c; // c>1
}

void graph_gnm_slow(graph_t g, unsigned long m);
extern long long vn_rng_calls;

int main(int argc, char* argv[]) {
  int i,j,k,mc,*cs,n=100,m,ntrials=1;
  double c,cmin=1.0,cmax=6.0,ex;
  graph_t g;
  if (argc>1) n=atoi(argv[1]);
  if (argc>2) cmin=atof(argv[2]);
  if (argc>3) cmax=atof(argv[3]);
  if (argc>4) ntrials=atoi(argv[4]);
  g=graph_new(n);
  for (k=0; k<100; k++) {
    c=cmin+(cmax-cmin)*k/100.0; // average degree
    ex=expected(n,c);
    m=(int)floor(0.5*c*n);
    for (i=0; i<ntrials; i++) {
      graph_gnm(g,m);
      //graph_gnm_slow(g,m);
      if (0) { // complete list of clusters
        cs=graph_cluster_sizes(g);
        j=0;
        while (cs[j]>1) { // clusters of 2 or more nodes
          printf("%d %d\n",j,cs[j]);
          j++;
        }
        free(cs);
      } else { // size of largest cluster
        mc=graph_max_cluster(g);
        printf("%g %d %g %g\n",c,mc,(mc/(double)n),ex);
      }
      fflush(stdout);
    }
  }
  printf("%lld RNG calls\n",vn_rng_calls);
  return 0;
}
