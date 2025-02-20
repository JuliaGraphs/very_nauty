// maximum degree of GRG(n,r=a*sqrt(log(n)/n))
// Penrose p 125
// for i in 1 2 3 4 5 6 7 8 9; do ./example_04 0.$i | grep nnodes | kw alpha lhs_on_rhs; done

#include <math.h>
#include "vn_graph.h"
#include "LambertW.c"

int main(int argc, char* argv[]) {
  int i,j,n;
  //double a=1.0,r,ns[]={100,150,200,250,300,500,1000,1500,2000};
  double a=1.0,r,ns[]={1000,1500,2000};
  double alpha,lhs,rhs;
  char chi_lbl[100];
  graph_t g;
  histogram_t h_Delta;
  if (argc>1) a=atof(argv[1]);
  alpha=a*a*M_PI; // FIXME *2?
  for (j=0; j<sizeof(ns)/sizeof(double); j++) {
    n=ns[j];
    g=graph_new(n);
    r=a*sqrt(log(n)/n);
    snprintf(chi_lbl,100,"Delta(grg(%d,%.3f))",n,r);
    h_Delta=histogram_new(chi_lbl);
    i=0;
    do {
      graph_grg_torus(g,r);
      histogram_add(h_Delta,graph_max_degree(g));
      if (n>0 && ++i%1000==0) {
        fprintf(stderr,"Not converged yet, i=%d\n",i);
        histogram_write(stderr,h_Delta);
      }
    } while(i<100 || !histogram_mean_done(h_Delta,1.0e-2));
    histogram_show(h_Delta);
    printf("n=%d r=%g mean_Delta_on_logn=%g ",n,r,histogram_mean(h_Delta)/log(n));
    printf("50_percentile=%g\n",histogram_quantile(h_Delta,0.5));
    lhs=histogram_mode(h_Delta)/(M_PI*n*r*r);
    //lhs=histogram_mean(h_Delta)/(M_PI*n*r*r);
    //lhs=histogram_quantile(h_Delta,0.5)/(M_PI*n*r*r);
    rhs=exp(1.0+LambertW((1.0/alpha-1.0)/M_E));
    printf("thm (6.14) a=%g alpha=%g nnodes=%d lhs=%g rhs=%g lhs_on_rhs=%g\n",a,alpha,n,lhs,rhs,lhs/rhs);
    fflush(stdout);
    histogram_clear(h_Delta);
    graph_clear(g);
  }
  return 0;
}

