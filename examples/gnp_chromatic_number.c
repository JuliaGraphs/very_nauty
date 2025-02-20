// KMB 2006 Aug 22
// make gnp_chromatic_number && time nice ./gnp_chromatic_number > gnp_chi.dat
// make gnp_chromatic_number && time nice ./gnp_chromatic_number | graph -Tx -C -g3 -x 0 16 -y 0 7 -m0 -S 16

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  int i,j,k=0,chi,to=0,n=200;
  double d,p;
  graph_t g;
  if (argc>1) n=atoi(argv[1]);
  g=graph_new(n);
  for (j=0; j<=100; j+=2) { // p loop
    d=0.16*j;
    p=d/n;
    //fprintf(stderr,"p=%g\n",p);
    for (i=0; i<10; i++) {
      k++;
      graph_gnp(g,p);
      chi=graph_chromatic_number(g,100);
      if (chi) printf("%g %d\n",d,chi); else to++;
    }
    fflush(stdout);
  }
  graph_clear(g);
  fprintf(stderr,"%d timeouts out of %d trials\n",to,k);
  return 0;
}

