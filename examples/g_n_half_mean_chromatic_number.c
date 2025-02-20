// KMB 2006 Nov 15
// make g_n_half_mean_chromatic_number && time nice ./g_n_half_mean_chromatic_number > gnp_chi.dat
// time nice ./g_n_half_mean_chromatic_number | graph -Tx -C -g3 -x 0 100 -y 0 20 -m1 -S 16 0.01

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  int i,k=0,chi,to,n,nmax=100;
  double a,p=0.5;
  graph_t g;
  if (argc>1) nmax=atoi(argv[1]);
  for (n=5; n<=nmax; n+=5) { // n loop
    g=graph_new(n);
    a=0.0;
    to=0;
    for (i=0; i<20; i++) {
      k++;
      graph_gnp(g,p);
      chi=graph_chromatic_number(g,1000);
      if (chi) a+=(chi-a)/(i+1); else to++;
    }
    if (a>0.0) printf("%d %g\n",n,a);
    fflush(stdout);
    graph_clear(g);
    fprintf(stderr,"n=%d: %d timeouts out of %d trials\n",n,to,k);
  }
  return 0;
}

