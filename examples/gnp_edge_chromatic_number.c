// KMB 2006 Nov 23
// distribution of edge chromatic number for G(n,d/n) for n fixed.
// make gnp_edge_chromatic_number && time nice ./gnp_edge_chromatic_number > gnp_chi_dash.dat
// time nice ./gnp_edge_chromatic_number | graph -Tx -C -g3 -x 0 16 -y 0 35 -m0 -S 16
// time nice ./gnp_edge_chromatic_number | graph -Tps -C -g3 -x 0 16 -y 0 35 -m0 -S 16 -L"edge chromatic number for G(n,d/n)" -Xd -Y"\*x'" >| gnp_edge_chromatic_number.ps
// mcopy -o -v gnp_edge_chromatic_number.c a:

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  int i,j,dj=2,k=0,chi,to=0,n=200,dmax=16;
  double d,p,mm,*meanmaxdeg;
  graph_t g;
  if (argc>1) n=atoi(argv[1]);
  meanmaxdeg=malloc(100*sizeof(double));
  g=graph_new(n);
  for (j=0; j<=100; j+=dj) { // d loop
    d=0.01*j*dmax;
    p=d/n;
    mm=0.0;
    for (i=0; i<50; i++) {
      k++;
      graph_gnp(g,p);
      mm+=(graph_max_degree(g)-mm)/(i+1.0);
      chi=graph_edge_chromatic_number(g,1000);
      if (chi) printf("%g %d\n",d,chi); else to++;
    }
    meanmaxdeg[j]=mm;
    fflush(stdout);
  }
  graph_clear(g);
  fprintf(stderr,"%d timeouts out of %d trials\n",to,k);
  printf("\n\n#m=2,S=4\n");
  for (j=0; j<=100; j+=dj) { // print mean Deltas
    d=0.01*j*dmax;
    printf("%g %g\n",d,meanmaxdeg[j]);
  }
  return 0;
}

