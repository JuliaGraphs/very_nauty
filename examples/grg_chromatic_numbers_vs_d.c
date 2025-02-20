// KMB 2006 Nov 27
// distribution of edge chromatic number for GRG(n,r) for n fixed, versus mean degree d.
// n=100; make grg_chromatic_numbers_vs_d && time nice ./grg_chromatic_numbers_vs_d $n | graph -Tx -C -g3 -x 0.2 0.8 -y 0 5 -S 0 0.01 -L"chromatic numbers for GRG($n,r)" -Xd -Y"\*x, \*x'"
// make grg_chromatic_numbers_vs_d && time nice ./grg_chromatic_numbers_vs_d | graph -Tps -C -g3 -m0 -S 16 0.01 -L"chromatics number for GRG(100,r)" -Xd -Y"\*x'" >| grg_chromatic_numbers_vs_d.ps
// mcopy -o -v grg_chromatic_numbers_vs_d.c a:

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  const int nj=100,timeout=100,nit=10000;
  int i,j,dj=1,k=0,chi,chi_dash,to=0,nnodes=100;
  double d,r,mm,*medchi,*medchid,dmin=0.2,dmax=0.8;
  FILE *foo=fopen("foo","w");
  graph_t g;
  histogram_t hchi,hchid;
  if (argc>1) nnodes=atoi(argv[1]);
  medchi =malloc(nj*sizeof(double));
  medchid=malloc(nj*sizeof(double));
  g=graph_new(nnodes);
  for (j=0; j<nj; j+=dj) { // d=mean degree loop
    hchi =histogram_new("chi(GRG)");
    hchid=histogram_new("chi_dash(GRG)");
    d=dmin+(dmax-dmin)*j/nj;
    r=sqrt(d/M_PI/(nnodes-1.0));
    mm=0.0;
    for (i=0; i<nit; i++) {
      k++;
      graph_grg_torus(g,r);
      mm+=(graph_max_degree(g)-mm)/(i+1.0);
      chi=graph_chromatic_number(g,timeout);
      if (chi>0) {
        printf("#m=0,S=1\n");
        printf("%g %d\n",d,chi); 
        histogram_add(hchi,chi);
      } 
        else to++;
      chi_dash=graph_edge_chromatic_number(g,timeout);
      if (chi_dash>=0) {
        printf("\n#m=1,S=2\n");
        printf("%g %d\n",d,chi_dash); 
        histogram_add(hchid,chi_dash);
      } 
        else to++;
    }
    medchi[j] =histogram_quantile(hchi,0.5);
    medchid[j]=histogram_quantile(hchid,0.5);
    if (medchid[j]>20) fprintf(stderr,"%g\n",medchid[j]), histogram_write(foo,hchid), exit(2);
    fflush(stdout);
    histogram_clear(hchi);
    histogram_clear(hchid);
  }
  graph_clear(g);
  printf("\n");
  fprintf(stderr,"%d timeouts out of %d trials\n",to,k);
  for (j=0; j<nj; j+=dj) {
    d=dmin+(dmax-dmin)*j/nj;
    printf("%g %g\n",d,medchi[j]);
  }
  printf("\n\n#m=3,S=0\n");
  for (j=0; j<nj; j+=dj) {
    d=dmin+(dmax-dmin)*j/nj;
    printf("%g %g\n",d,medchid[j]);
  }
  return 0;
}

