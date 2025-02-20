// KMB 2006 Nov 25
// distribution of edge chromatic number for GRG(n,r) for n fixed.
// make grg_edge_chromatic_number && time nice ./grg_edge_chromatic_number | graph -Tx -C -g3 -x 0 0.2 -y 0 25 -m0 -S 16 0.01 -L"edge chromatic number for GRG(n,r)" -Xr -Y"\*x'"
// make grg_edge_chromatic_number && time nice ./grg_edge_chromatic_number | graph -Tps -C -g3 -m0 -S 16  0.01 -L"edge chromatic number for GRG(100,r)" -Xr -Y"\*x'" >| grg_edge_chromatic_number.ps
// mcopy -o -v grg_edge_chromatic_number.c a:

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  int i,j,dj=2,k=0,chi_dash,to=0,n=100;
  double r,mm,*lq,*uq,rmax=0.2;
  graph_t g;
  histogram_t h;
  if (argc>1) n=atoi(argv[1]);
  lq=malloc(100*sizeof(double));
  uq=malloc(100*sizeof(double));
  g=graph_new(n);
  for (j=0; j<=100; j+=dj) { // r loop
    h=histogram_new("chi_dash(GRG)");
    r=rmax*j/100.0;
    mm=0.0;
    for (i=0; i<200; i++) {
      k++;
      graph_grg_torus(g,r);
      //graph_show(g);
      mm+=(graph_max_degree(g)-mm)/(i+1.0);
      chi_dash=graph_edge_chromatic_number(g,10);
      if (chi_dash>=0) {
        printf("%g %d\n",r,chi_dash); 
        histogram_add(h,chi_dash);
      } 
        else to++;
    }
    lq[j]=histogram_quantile(h,0.1);
    uq[j]=histogram_quantile(h,0.9);
    fflush(stdout);
    histogram_clear(h);
  }
  graph_clear(g);
  fprintf(stderr,"%d timeouts out of %d trials\n",to,k);
  printf("\n\n#m=2,S=0\n");
  for (j=0; j<=100; j+=dj) {
    r=0.01*j*rmax;
    printf("%g %g\n",r,lq[j]);
  }
  printf("\n\n#m=3,S=0\n");
  for (j=0; j<=100; j+=dj) {
    r=0.01*j*rmax;
    printf("%g %g\n",r,uq[j]);
  }
  return 0;
}

