// KMB 2007 Jan 23
// 

#include <math.h>
#include "vn_graph.h"

typedef struct {
  histogram_t chi;
} graph_data_t;

int f(graph_t g, int i, void *cd) {
  // timeout will return -1, which histogram_add ignores...
  histogram_add(((graph_data_t*)cd)->chi,graph_chromatic_number(g,100)); 
  if (i>=99) {
    return 0; // done
  }
  return 1; // continue
}

int main(int argc, char* argv[]) {
  unsigned int i,n=400,nd=10;
  double d,p;
  graph_data_t cd;
  FILE* fchi=fopen("gnp_chi_01.his","w");
  cd.chi=histogram_new("chromatic number of G(n,d/n)");
  if (argc>1) n=atoi(argv[1]);  // n
  for (i=0; i<nd; i++) {
    d=(16.0*i)/nd;
    d=4.0+i/(double)nd; // zoom 4..5
    p=d/n;
    fprintf(fchi,"d=%g\n",d);
    histogram_empty(cd.chi);
    graph_gnp_iterator(n,p,f,(void*)&cd);
    histogram_write(fchi,cd.chi);
    fflush(fchi);
  }
  histogram_clear(cd.chi);
  fclose(fchi);
  return 0;
}

