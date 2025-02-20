//

#include <math.h>
#include "vn_graph.h"

int main(int argc, char* argv[]) {
  int i,j,n,chi;
  //double a=1.0,r,ns[]={5,10,12,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,180,200,250,300};
  double a=1.0,r,ns[]={100};
  char chi_lbl[100],ub_lbl[100];
  graph_t g;
  histogram_t h_chi,h_ub;
  if (argc>1) a=atof(argv[1])/100.0;
  for (j=0; j<sizeof(ns)/sizeof(double); j++) {
    n=ns[j];
    g=graph_new(n);
    r=a*sqrt(log(n)/n);
    snprintf(chi_lbl,100,"chi(grg(%d,%.3f))=",n,r);
    h_chi=histogram_new(chi_lbl);
    snprintf(ub_lbl,100,"chi(grg(%d,%.3f))<=",n,r);
    h_ub=histogram_new(ub_lbl);
    i=0;
    do {
      graph_grg_torus(g,r);
      chi=graph_chromatic_number(g,20);
      if (chi>0)
        histogram_add(h_chi,chi);
      else {
        fprintf(stderr,"timeout i=%d\n",i);
        histogram_add(h_ub,graph_sequential_color_repeat(g,1000));
      }
      if (n>0 && ++i%100==0) {
        fprintf(stderr,"Not converged yet, i=%d\n",i);
        histogram_write(stderr,h_chi);
        histogram_write(stderr,h_ub);
      }
    } while(i<50 || !histogram_mean_done(h_chi,0.5e-2));
    histogram_show(h_chi);
    histogram_show(h_ub);
    printf("a=%g nnodes=%d r=%g mean_chi_on_logn=%g ",a,n,r,histogram_mean(h_chi)/log(n));
    printf("50_percentile=%g\n",histogram_quantile(h_chi,0.5));
    fflush(stdout);
    histogram_clear(h_chi);
    histogram_clear(h_ub);
    graph_clear(g);
  }
  return 0;
}

