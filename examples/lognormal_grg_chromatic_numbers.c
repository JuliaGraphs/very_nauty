// KMB 2007 Jan 10
// distribution of edge chromatic number for ln_GRG(n,r) for n fixed.
// make lognormal_grg_chromatic_numbers && time nice ./lognormal_grg_chromatic_numbers | graph -Tx -C -g3 -x 0 0.06 -y 0 5 -S 0 0.01 -L"chromatic numbers for ln_GRG(100,r,alpha)" -Xr -Y"\*x, \*x'"
// mcopy -o -v lognormal_grg_chromatic_numbers.c a:

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  int i,j,dj=2,k=0,chi,chi_dash,to=0,n=100;
  double r,mm,*medchi,*medchid,rmax=0.06,alpha=5.0e-1;
  graph_t g;
  histogram_t hchi,hchid;
  if (argc>1) n=atoi(argv[1]);
  medchi =malloc(100*sizeof(double));
  medchid=malloc(100*sizeof(double));
  g=graph_new(n);
  for (j=0; j<=100; j+=dj) { // r loop
    hchi =histogram_new("chi(ln_GRG)");
    hchid=histogram_new("chi_dash(ln_GRG)");
    r=rmax*j/100.0;
    mm=0.0;
    for (i=0; i<100; i++) {
      k++;
      graph_lognormal_grg_torus(g,r,alpha);
      mm+=(graph_max_degree(g)-mm)/(i+1.0);
      chi=graph_chromatic_number(g,100);
      if (chi>0) {
        printf("#m=0,S=1\n");
        printf("%g %d\n",r,chi); 
        histogram_add(hchi,chi);
      } 
        else to++;
      chi_dash=graph_edge_chromatic_number(g,100);
      if (chi_dash>=0) {
        printf("\n#m=1,S=2\n");
        printf("%g %d\n",r,chi_dash); 
        histogram_add(hchid,chi_dash);
      } 
        else to++;
    }
    medchi[j] =histogram_quantile(hchi,0.5);
    medchid[j]=histogram_quantile(hchid,0.5);
    fflush(stdout);
    histogram_clear(hchi);
    histogram_clear(hchid);
  }
  graph_clear(g);
  fprintf(stderr,"%d timeouts out of %d trials\n",to,k);
  printf("\n\n#m=2,S=0\n");
  for (j=0; j<=100; j+=dj) {
    r=0.01*j*rmax;
    printf("%g %g\n",r,medchi[j]);
  }
  printf("\n\n#m=3,S=0\n");
  for (j=0; j<=100; j+=dj) {
    r=0.01*j*rmax;
    printf("%g %g\n",r,medchid[j]);
  }
  return 0;
}

