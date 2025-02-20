// KMB 2006 Nov 15
// make gnp_mean_chromatic_number && time nice ./gnp_mean_chromatic_number > gnp_chi.dat
// time nice ./gnp_mean_chromatic_number | graph -Tx -C -g3 -x 0 100 -y 0 20 -m1 -S 16 0.01

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  int j,k=0,chi,to,n,nmax=100;
  double a,aa,p=0.5,ss,eps=5e-3;
  graph_t g;
  if (argc>1) nmax=atoi(argv[1]);
  if (argc>2) p=atof(argv[2]);
  for (n=5; n<=nmax; n+=5) { // n loop
    g=graph_new(n);
    ss=a=0.0;
    to=0;
    j=k=0; 
    do {
      k++;
      graph_gnp(g,p);
      chi=graph_chromatic_number(g,1000);
      if (chi) {
        j++;
        aa=a;
        a+=(chi-a)/j; 
        ss+=(chi-a)*(chi-aa);
        if (j>10 && sqrt(ss/j)/sqrt(j)<eps*a) goto done; // stderr
      }
        else to++;
    } while ((n<20 && k<10000) || (n<50 && k<1000) || (n>=50 && k<100));
    fprintf(stderr,"n=%d: mean not converged, k=%d, mean=%g reldev=%g\n",n,k,a,sqrt(ss/j)/a);
    done:
    if (a>0.0) printf("%d %g\n",n,a);
    fflush(stdout);
    graph_clear(g);
    if (to>0) fprintf(stderr,"n=%d: %d timeouts out of %d trials\n",n,to,k);
  }
  return 0;
}

